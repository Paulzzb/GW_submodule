# ISDF inline debug controls (`+isdf/+debug`)

## Session snapshot (recommended)

Debug flags are read from **`config.ISDF` once** and cached for the MATLAB session:

- **`isdftest.debug.init_from_config(config)`** is called from **`service_driver`** at startup, and from **`isdftest.adaptive.run_adaptiveisdf`** after relay restore (paths that skip `service_driver`).
- After init, call sites use **`isdftest.debug.on(tag)`** and **`isdftest.debug.react(...)`** with **no `config` argument**.
- **`isdftest.debug.clear()`** resets the cache; it is registered in **`service_reset_persistent`** alongside other ISDF resets.

| Field (`config.ISDF`) | Default | Meaning |
|------------------------|---------|--------|
| `debug_checks` | `false` | Master switch for tagged inline checks. |
| `debug_level` | `'error'` | `'error'` → `error()` on violation; `'warn'` / `'warning'` → `warning()`. |
| `debug_tags` | `[]` | Empty → all tags when `debug_checks` is true; nonempty → only listed tags (case-insensitive). |

Defaults are set in `GW/input/default_param_values.m`.

---

## API

### `isdftest.debug.init_from_config(config)`

Snapshot `debug_*` from `config.ISDF`. Safe if `ISDF` is missing (checks stay off).

### `isdftest.debug.clear()`

Clear the snapshot (checks off until the next `init_from_config`).

### `tf = isdftest.debug.on(tag)`

Use **outside** expensive work (`if isdftest.debug.on(tag) ... end`). Uses **only** the session cache. If `init_from_config` was never called, returns `false`.

### `isdftest.debug.react(violation, msg [, tagForId])`

If `violation` is true, emit `warning` or `error` per cached `debug_level`. Optional `tagForId` tail for the warning identifier (sanitized to `[a-zA-Z0-9_]`).

### Overloads (optional, for scripts without `service_driver`)

- **`isdftest.debug.on(config, tag)`** — read `config.ISDF` directly; does **not** update the cache.
- **`isdftest.debug.level(config)`** / **`isdftest.debug.level()`** — level from config or cache.
- **`isdftest.debug.react(config, violation, msg [, tagForId])`**
- **`isdftest.debug.check(tag, fh, msg [, id])`** or **`isdftest.debug.check(config, tag, fh, msg [, id])`**

### Internal

- **`isdftest.debug.cache`** — persistent store; not intended for general use.

---

## Usage (Mode B + react)

```matlab
if isdftest.debug.on('rsymm/bundle_refresh')
  % expensive checks only when enabled
  isdftest.debug.react(norm(a - b) > tol, 'message', 'my_check_id');
end
```

### Tags in use

| Tag | Location |
|-----|----------|
| `coeff/coeff_coarse_wf_extract` | `coeff_coarse_wf_extract.m` — coincident coarse vs fine WF |
| `rsymm/bundle_refresh` | `bundle_refresh.m` — one gate for: torus rotation check, discrete `R_rot` map on bundle rows, full-band WF bundle/sampling checks, `ib=3:5` sampling-vs-rotated-path check, and **post-`save2mod`** `bs_new` tests (`ib=3:5`, `S_q` fine/coarse vs direct) (all skipped when debug off) |

---

## Scripts without `service_driver`

Call **`isdftest.debug.init_from_config(config)`** after you have the same `config` struct you would pass to `service_driver` (e.g. `load(..., ''config'')`).

---

**Path:** `GW/service/+isdf/+debug/ISDF_debug.md`  
**Chinese:** `GW/service/+isdf/+debug/ISDF_debug_zh.md`  
**Last updated:** 2026-05-02
