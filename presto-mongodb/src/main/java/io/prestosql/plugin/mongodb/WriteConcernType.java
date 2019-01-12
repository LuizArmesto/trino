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
package io.prestosql.plugin.mongodb;

import com.mongodb.WriteConcern;

import static java.util.Objects.requireNonNull;

public enum WriteConcernType
{
    ACKNOWLEDGED(WriteConcern.ACKNOWLEDGED),
    FSYNC_SAFE(WriteConcern.FSYNC_SAFE),
    FSYNCED(WriteConcern.FSYNCED),
    JOURNAL_SAFEY(WriteConcern.JOURNAL_SAFE),
    JOURNALED(WriteConcern.JOURNALED),
    MAJORITY(WriteConcern.MAJORITY),
    NORMAL(WriteConcern.NORMAL),
    REPLICA_ACKNOWLEDGED(WriteConcern.REPLICA_ACKNOWLEDGED),
    REPLICAS_SAFE(WriteConcern.REPLICAS_SAFE),
    SAFE(WriteConcern.SAFE),
    UNACKNOWLEDGED(WriteConcern.UNACKNOWLEDGED);

    private final WriteConcern writeConcern;

    WriteConcernType(WriteConcern writeConcern)
    {
        this.writeConcern = requireNonNull(writeConcern, "writeConcern is null");
    }

    public WriteConcern getWriteConcern()
    {
        return writeConcern;
    }
}
