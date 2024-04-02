/*
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package io.trino.plugin.geospatial;

import com.esri.core.geometry.Point;
import com.esri.core.geometry.ogc.OGCGeometry;
import com.esri.core.geometry.ogc.OGCGeometryCollection;
import com.esri.core.geometry.ogc.OGCLineString;
import com.esri.core.geometry.ogc.OGCPoint;
import com.esri.core.geometry.ogc.OGCPolygon;
import com.google.errorprone.annotations.FormatMethod;
import com.uber.h3core.H3Core;
import com.uber.h3core.exceptions.H3Exception;
import com.uber.h3core.util.LatLng;
import io.airlift.slice.Slice;
import io.trino.geospatial.GeometryType;
import io.trino.spi.TrinoException;
import io.trino.spi.block.Block;
import io.trino.spi.block.BlockBuilder;
import io.trino.spi.function.Description;
import io.trino.spi.function.ScalarFunction;
import io.trino.spi.function.SqlNullable;
import io.trino.spi.function.SqlType;
import io.trino.spi.type.StandardTypes;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static io.trino.geospatial.serde.GeometrySerde.deserialize;
import static io.trino.plugin.geospatial.GeometryType.GEOMETRY_TYPE_NAME;
import static io.trino.spi.StandardErrorCode.INVALID_FUNCTION_ARGUMENT;
import static io.trino.spi.type.BigintType.BIGINT;
import static java.lang.String.format;

/**
 * A set of functions to convert between geometries and H3 cells.
 *
 * @see <a href="https://h3geo.org/docs">https://h3geo.org/docs</a>
 * for the description of the H3 geospatial indexing system.
 */
public final class H3Functions
{
    public static final int MAX_RESOLUTION_LEVEL = 15;
    static final String TYPE_ARRAY_BIGINT = "ARRAY(BIGINT)";
    private static final String RESOLUTION_LEVEL_TOO_SMALL = "Resolution level must be >= 0";
    private static final String RESOLUTION_LEVEL_TOO_LARGE = "Resolution level must be <= " + MAX_RESOLUTION_LEVEL;

    public static H3Core h3;
    public static Map<Integer, Double> cellSizeMap;

    static {
        try {
            h3 = H3Core.newInstance();
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
        // init cellSizeMap
        cellSizeMap = new HashMap<>();
        cellSizeMap.put(0, 1257.1651279870703);
        cellSizeMap.put(1, 415.1068893923832);
        cellSizeMap.put(2, 147.84174424387626);
        cellSizeMap.put(3, 54.530253595150185);
        cellSizeMap.put(4, 20.426305676421265);
        cellSizeMap.put(5, 7.693156379723934);
        cellSizeMap.put(6, 2.9039922554770365);
        cellSizeMap.put(7, 1.097050477001042);
        cellSizeMap.put(8, 0.4145696618254594);
        cellSizeMap.put(9, 0.15668127143434735);
        cellSizeMap.put(10, 0.059218394287200805);
        cellSizeMap.put(11, 0.022382219491694798);
        cellSizeMap.put(12, 0.008459648046313144);
        cellSizeMap.put(13, 0.0031974550140400557);
        cellSizeMap.put(14, 0.0012085926447938334);
        cellSizeMap.put(15, 0.0004567448929842399);
    }

    private H3Functions() {}

    public static double latLngDistance(double lat1, double lat2, double lon1, double lon2)
    {
        final int radius = 6371; // Radius of the earth

        double latDistance = Math.toRadians(lat2 - lat1);
        double lonDistance = Math.toRadians(lon2 - lon1);
        double a = Math.sin(latDistance / 2) * Math.sin(latDistance / 2)
                + Math.cos(Math.toRadians(lat1)) * Math.cos(Math.toRadians(lat2))
                * Math.sin(lonDistance / 2) * Math.sin(lonDistance / 2);
        double c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
        return radius * c;
    }

    public static List<Long> approxPathCells(OGCPoint sc, OGCPoint ec, int level, boolean safe)
    {
        double approxDistance = latLngDistance(sc.Y(), ec.Y(), sc.X(), ec.X());
        int itCellNum = Math.max((int) (approxDistance / cellSizeMap.get(level) * 3), 1);
        double xUnit = (ec.X() - sc.X()) / itCellNum;
        double yUnit = (ec.Y() - sc.Y()) / itCellNum;
        // keep order of cells
        Set<Long> cells = new HashSet<>();
        for (int it = 0; it <= itCellNum; ++it) {
            cells.addAll(pointToH3Cells(new OGCPoint(new Point(sc.X() + xUnit * it, sc.Y() + yUnit * it), sc.esriSR), level, false));
        }
        // truncate the duplicate case where 1 cell connect to > 1 neighbor cell.
        List<Long> approxCells = new ArrayList<>(cells);
        // clear cells, put the first cell as start
        cells.clear();
        cells.add(approxCells.get(0));
        Set<Long> closetNeighbors = new HashSet<>(h3.gridDisk(approxCells.get(0), 1));
        List<Long> candidates = new ArrayList<>();
        int it = 1;
        while (it < approxCells.size()) {
            // if it's one of the neighbors of the closet cell,
            long curCell = approxCells.get(it);
            if (closetNeighbors.contains(curCell)) {
                // do 1 step backward check
                candidates.add(curCell);
                it++;
            }
            else {
                // the curCell is not the neighbor of closest cell. means check for candidates ended.
                long pickedCandidate = 0;
                for (Long candidate : candidates) {
                    List<Long> candidateNeighbors = h3.gridDisk(candidate, 1);
                    if (candidateNeighbors.contains(curCell)) {
                        // find the cell that connect to both cur and closest, add into final list
                        pickedCandidate = candidate;
                        // if the candidate connect to both cur and next cell, means it's better to pick it
                        if (it + 1 < approxCells.size() && candidateNeighbors.contains(approxCells.get(it + 1))) {
                            break;
                        }
                    }
                }
                if (pickedCandidate != 0) {
                    cells.add(pickedCandidate);
                    candidates.clear();
                    closetNeighbors = new HashSet<>(h3.gridDisk(pickedCandidate, 1));
                }
                else if (!safe) {
                    // if not need to guarantee safe, randomly pick 1 candidate
                    cells.add(candidates.get(0));
                    cells.add(curCell);
                    closetNeighbors = new HashSet<>(h3.gridDisk(curCell, 1));
                    candidates.clear();
                    // if not raise error, make sure it always rotate
                    it++;
                }
                else {
                    // our estimation shouldn't generate two adjacent cells in list that are not neighbors.
                    // because we split the list by the min cell edge length / 3, and the least possible distance between the edges of 2 hexagon
                    // with 1 hexagon in between equals to min cell edge length. We'll throw error here if the edge case found and use it for corner case study.
                    // no matter what, move the iteration
                    throw new RuntimeException(
                        String.format(
                            "fail to approximate path between cell %d to %d, generated gap cell %d",
                            pointToH3Cells(sc, level, false).get(0),
                            pointToH3Cells(ec, level, false).get(0),
                            curCell));
                }
            }
        }
        if (!candidates.isEmpty()) {
            cells.add(candidates.get(0));
        }
        return new ArrayList<>(cells);
    }

    private static List<Long> polygonToH3Cells(OGCPolygon polygon, int resolution, boolean fullCover)
    {
        List<LatLng> shell = new ArrayList<>();
        List<List<LatLng>> holes = new ArrayList<>();

        OGCLineString exteriorRing = polygon.exteriorRing();
        for (int i = 0; i < exteriorRing.numPoints(); i++) {
            OGCPoint point = exteriorRing.pointN(i);
            shell.add(new LatLng(point.Y(), point.X()));
        }

        for (int i = 0; i < polygon.numInteriorRing(); i++) {
            OGCLineString interiorRing = polygon.interiorRingN(i);
            List<LatLng> hole = new ArrayList<>();
            for (int ii = 0; ii < exteriorRing.numPoints(); ii++) {
                OGCPoint point = interiorRing.pointN(ii);
                hole.add(new LatLng(point.Y(), point.X()));
            }
            if (!hole.isEmpty()) {
                holes.add(hole);
            }
        }
        // H3 polyfill only include hexagons with centroid within the polygon, we fix by generating the neighbors for the cells.
        List<Long> internalCells = h3.polygonToCells(shell, holes.isEmpty() ? null : holes, resolution);

        // if the internal cell is empty, means no cell with centroid within, we'll use centroid to add a seed
        if (internalCells.isEmpty()) {
            OGCPoint centroid = (OGCPoint) polygon.centroid();
            internalCells.add(h3.latLngToCell(centroid.Y(), centroid.X(), resolution));
        }

        List<Long> cells;
        if (fullCover) {
            cells = new ArrayList<>();
            // TODO
        }
        else {
            cells = internalCells;
        }

        return cells;
    }

    private static List<Long> linestringToH3Cells(OGCLineString linestring, int resolution, boolean fullCover)
    {
        Set<Long> cells = new HashSet<>();

        List<OGCPoint> points = new ArrayList<>();
        for (int i = 0; i < linestring.numPoints(); i++) {
            points.add(linestring.pointN(i));
        }

        // for each segment, generate coverage separately
        for (int pairE = 1; pairE < points.size(); pairE++) {
            OGCPoint cs = points.get(pairE - 1);
            OGCPoint ce = points.get(pairE);

            long cellS = pointToH3Cells(cs, resolution, fullCover).get(0);
            long cellE = pointToH3Cells(ce, resolution, fullCover).get(0);

            List<Long> pathCells;
            try {
                try {
                    pathCells = h3.gridPathCells(cellS, cellE);
                }
                catch (H3Exception e) {
                    pathCells = approxPathCells(cs, ce, resolution, true);
                }

                if (pathCells.size() <= 2 || !fullCover) {
                    cells.addAll(pathCells);
                }
                else {
                    // the shortest path from cs to ce not guarantee cover the line, for the missed part, we generate neighbors for the cells to cover.
                    // We exclude the cells that only adjacent to 1 cell in the path to remove some redundancy, but this won't remove all redundancy.
                    Set<Long> prevNeighbors = new HashSet<>(h3.gridDisk(pathCells.get(0), 1));
                    for (int it = 1; it < pathCells.size(); it++) {
                        Set<Long> neighbors = new HashSet<>(h3.gridDisk(pathCells.get(it), 1));
                        prevNeighbors.retainAll(neighbors);
                        cells.addAll(prevNeighbors);
                        prevNeighbors = neighbors;
                    }
                }
            }
            catch (Exception e) {
                if (fullCover) {
                    // TODO
                }
                else {
                    // if not full cover, we'll throw error, cause for getting the shortest path, we suppose should have no error
                    throw e;
                }
            }
        }

        return cells.stream().toList();
    }

    private static List<Long> pointToH3Cells(OGCPoint point, int resolution, boolean fullCover)
    {
        return List.of(h3.latLngToCell(point.Y(), point.X(), resolution));
    }

    @Description("Cover the geometry by H3 cell IDs with the given resolution.")
    @ScalarFunction("h3_cell_ids")
    @SqlNullable
    @SqlType(H3Functions.TYPE_ARRAY_BIGINT)
    public static Block h3CellIds(@SqlType(GEOMETRY_TYPE_NAME) Slice polygonSlice, @SqlType(StandardTypes.INTEGER) long resolution, @SqlType(StandardTypes.BOOLEAN) boolean fullCover)
    {
        checkResolutionLevel(resolution);

        OGCGeometry geometry = deserialize(polygonSlice);
        if (geometry.isEmpty()) {
            return null;
        }

        HashSet<Long> cellIds = new HashSet<>();

        List<OGCGeometry> geoms = new ArrayList<>();
        if (geometry instanceof OGCGeometryCollection collection) {
            for (int i = 0; i < collection.numGeometries(); i++) {
                geoms.add(collection.geometryN(i));
            }
        }

        // first flat all GeometryCollection implementations into single Geometry types (Polygon, LineString, Point)
        for (OGCGeometry geom : geoms) {
            GeometryType type = GeometryType.getForEsriGeometryType(geom.geometryType());
            if (type == GeometryType.POLYGON) {
                cellIds.addAll(polygonToH3Cells((OGCPolygon) geom, (int) resolution, fullCover));
            }
            else if (type == GeometryType.LINE_STRING) {
                cellIds.addAll(linestringToH3Cells((OGCLineString) geom, (int) resolution, fullCover));
            }
            else if (type == GeometryType.POINT) {
                cellIds.addAll(pointToH3Cells((OGCPoint) geom, (int) resolution, fullCover));
            }
            else {
                // if not type of polygon, point or lienSting, we cover its MBR
                cellIds.addAll(polygonToH3Cells((OGCPolygon) geom.envelope(), (int) resolution, fullCover));
            }
        }

        BlockBuilder blockBuilder = BIGINT.createFixedSizeBlockBuilder(cellIds.size());
        for (Long cell : cellIds) {
            BIGINT.writeLong(blockBuilder, cell);
        }
        return blockBuilder.build();
    }

    private static void checkResolutionLevel(long resolution)
    {
        checkCondition(resolution > 0, RESOLUTION_LEVEL_TOO_SMALL);
        checkCondition(resolution <= MAX_RESOLUTION_LEVEL, RESOLUTION_LEVEL_TOO_LARGE);
    }

    private static void checkCondition(boolean condition, String message)
    {
        checkCondition(condition, "%s", message);
    }

    @FormatMethod
    private static void checkCondition(boolean condition, String formatString, Object... args)
    {
        if (!condition) {
            throw new TrinoException(INVALID_FUNCTION_ARGUMENT, format(formatString, args));
        }
    }
}
