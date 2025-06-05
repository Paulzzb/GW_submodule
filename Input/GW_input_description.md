<!-- <pre lang="markdown"> &lt;style&gt; a { color: #0077cc} a:hover { text-decoration: underline; } &lt;/style&gt; </pre> -->

# GWOptions Input File Description

This document describes all supported input blocks and parameters used to configure the GWOptions framework. Parameters are grouped by block (namelist-style) and include descriptions, expected types, and default values.

---

## 🔎 Alphabetical Parameter Index

### &ISDF

[`isisdf`](#isisdf) |
[`isdf_ratio`](#isdf-ratio) |
[`isdf_ratio_type1`](#isdf-ratio-type1) |
[`isdf_ratio_type2`](#isdf-ratio-type2) |
[`isdf_ratio_type3`](#isdf-ratio-type3) |

---

## &ISDF Block

Controls the behavior of the ISDF (Interpolative Separable Density Fitting) approximation used in the GW workflow.

| Parameter            | Type    | Default      | Description                                                   |
|----------------------|---------|--------------|---------------------------------------------------------------|
| <a name="isdf-ratio"></a>`isdf_ratio`         | float   | `1.0`        | Global ratio applied to all ISDF-related components           |
| <a name="isdf-ratio-type1"></a>`isdf_ratio_type1`   | float   | `isdf_ratio` | Override ratio for ISDF type 1 operations                    |
| <a name="isdf-ratio-type2"></a>`isdf_ratio_type2`   | float   | `isdf_ratio` | Override ratio for ISDF type 2 operations                    |
| <a name="isdf-ratio-type3"></a>`isdf_ratio_type3`   | float   | `isdf_ratio` | Override ratio for ISDF type 3 operations                    |

> If any of the `isdf_ratio_typeX` parameters are omitted, they inherit the value of `isdf_ratio`.

---

## &FREQUENCY Block

_Future content to be added here._

---

## &CUTOFFS Block

_Future content to be added here._

---

## &CONTROL Block

_Future content to be added here._

---

## &SYSTEM Block

_Future content to be added here._

---

## Notes

- All parameter names are case-insensitive.
- Input follows a namelist-style format (`&BLOCK ... /`).
- Comments starting with `!` or `#` are ignored.
