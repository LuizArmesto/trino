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
package io.prestosql.plugin.tpcds;

import com.facebook.presto.spi.ConnectorTableLayoutHandle;
import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;

import static java.util.Objects.requireNonNull;

public class TpcdsTableLayoutHandle
        implements ConnectorTableLayoutHandle
{
    private final TpcdsTableHandle table;

    @JsonCreator
    public TpcdsTableLayoutHandle(@JsonProperty("table") TpcdsTableHandle table)
    {
        this.table = requireNonNull(table, "table is null");
    }

    @JsonProperty
    public TpcdsTableHandle getTable()
    {
        return table;
    }

    @Override
    public String toString()
    {
        return table.getTableName();
    }
}
