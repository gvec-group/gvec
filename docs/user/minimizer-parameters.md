:::{dropdown} `whichInitEquilibrium` (mandatory)
:name: whichInitEquilibrium

How initial guess is computed.
Either from boundary and axis parameters, or from VMEC file.

**Required Input**

**Type:** `integer`

**Allowed Values:**
:::{list-table}
*   - `0`
    - from axis and boundary parameters
*   - `1`
    - from VMEC file, **needs:** [`VMECwoutfile`](<project:#VMECwoutfile>),[`VMECwoutfile_format`](<project:#VMECwoutfile_format>)
:::

:::
:::{dropdown} `VMECwoutfile`
:name: VMECwoutfile

full file name of vmec solution file, either as netcdf or as nemec output,
see [`VMECwoutfile_format`](<project:#VMECwoutfile_format>)

**Required if:** [`whichInitEquilibrium=1`](<project:#whichInitEquilibrium>)

**Type:** `string`

:::
