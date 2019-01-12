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
package io.prestosql.type;

import com.facebook.presto.spi.block.Block;
import com.facebook.presto.spi.block.BlockBuilder;

import static com.facebook.presto.spi.type.IntegerType.INTEGER;

public class TestIntegerType
        extends AbstractTestType
{
    public TestIntegerType()
    {
        super(INTEGER, Integer.class, createTestBlock());
    }

    public static Block createTestBlock()
    {
        BlockBuilder blockBuilder = INTEGER.createBlockBuilder(null, 15);
        INTEGER.writeLong(blockBuilder, 1111);
        INTEGER.writeLong(blockBuilder, 1111);
        INTEGER.writeLong(blockBuilder, 1111);
        INTEGER.writeLong(blockBuilder, 2222);
        INTEGER.writeLong(blockBuilder, 2222);
        INTEGER.writeLong(blockBuilder, 2222);
        INTEGER.writeLong(blockBuilder, 2222);
        INTEGER.writeLong(blockBuilder, 2222);
        INTEGER.writeLong(blockBuilder, 3333);
        INTEGER.writeLong(blockBuilder, 3333);
        INTEGER.writeLong(blockBuilder, 4444);
        return blockBuilder.build();
    }

    @Override
    protected Object getGreaterValue(Object value)
    {
        return ((Long) value) + 1;
    }
}
