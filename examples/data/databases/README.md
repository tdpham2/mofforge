# Fictional database fixtures

Every identifier, DOI, CCDC number, pore dimension, stability value, and database
record in this folder is invented for teaching. These tables are not CSD exports
or measured/computed CoRE MOF data. No proprietary database records are bundled.

The tables use the column names consumed by mofforge. Unused columns can be
absent; the existing parsers supply empty fields. Four CSD records and four CoRE
records demonstrate the following relationships:

| CSD identifier | CoRE entry | Geometry copied at runtime |
| --- | --- | --- |
| DEMOZN | 2099[Zn][pcu]3[ASR]1 | Bundled IRMOF-1.cif |
| DEMOZN | 2099[Zn][pcu]3[FSR]1 | The **same** IRMOF-1.cif |
| DEMOZR | 2099[Zr][fcu]3[ASR]1 | Bundled UiO-66.cif; numerical metadata intentionally missing |
| DEMOMIS | 2099[Zn][pcu]3[ASR]2 | Intentionally absent |
| DEMONONE | No CoRE counterpart | None |

The ASR/FSR demonstration copies do not model different experimental activation
states. Their different metadata values exist solely to demonstrate filtering.
[structure_map.json](structure_map.json) records the mapping. Lessons stage
copies beneath their output directory and use the real database classes.

Change the lesson's explicit data-path arguments to work with your own datasets.
The [database guide](../../databases/README.md) explains missing values, multiple
processing variants, local structure resolution, and required real-data paths.
