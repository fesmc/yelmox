// Forcing builder: the forcing of one run as the equation of
// docs/_forcing-terms.md, built in the order of the equation. First, across
// the page, the explainer (what a snapshot, a term and a coefficient are, and
// an example), then the presets and the base par file. Left, top to
// bottom: the forcing time, the equation (the reference, the terms, then the
// coefficients they name; each with where its data comes from: the record of
// a coefficient, the fields of each snapshot), what comes after (transforms,
// orbit, sea level, prescribed ice). Right: the equation chain evaluated at a
// chosen time, the coefficients over the run, and the run written out (the
// forcing groups as a par file; the -p list against a base par file with the
// runme command).
//
// Words: the scalar α(τ) is a coefficient; in the namelist it is a group
// [coef_<name>], listed in [forcing] coefs and named by the key coef of a
// term (the three constants below).
//
// State: the forcing groups of the run, materialised (S.g). A group is
// created from its source on first use (the base par file, the library
// input/forcing_db.nml, or the schema template of input/yelmox_defaults.nml
// for a new name) and edited in place; the run is written out by comparing
// each group with its source. A preset or a base par file plus -p is read
// into the groups (import); the -p list is the difference (export).
// Data: forcing.json (scripts/forcing-builder/export.jl). The coefficient formulas
// mirror fesm-utils tsgen (TsGen.jl).

import * as Plot from "https://cdn.jsdelivr.net/npm/@observablehq/plot@0.6.16/+esm";
import katex from "https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.mjs";

const root = document.getElementById("fb");
let D = null;                       // forcing.json
let TEMPLATES = [];                 // [regex, template] (coef_, term_, snapshot_)
const S = { base: null, g: new Map(), extra: [], preset: null, presetKey: "", tsel: null, tab: "par", open: new Set(), runopts: "-r", outdir: "" };

// The namelist names of the coefficients: the group prefix, the list in
// [forcing], the key of a term. A rename of the namelist changes these three.
const COEF_PREFIX = "coef_", COEF_LIST_KEY = "coefs", TERM_COEF_KEY = "coef";
// the template of the coefficient groups in input/yelmox_defaults.nml
const COEF_T = COEF_PREFIX.replace(/_$/, "");

const FIELDS = ["tas", "pr", "to", "so"];
const SNAPFIELDS = ["tas", "pr", "zs", "to", "so"];
const FSYM = { tas: "T^{\\mathrm{sl}}", pr: "P^{\\mathrm{sl}}", to: "\\Theta", so: "S" };
const GAIN = { tas: "f_ta", to: "f_to", so: "f_so" };
const GREEK = { alpha: "\\alpha", beta: "\\beta", gamma: "\\gamma", delta: "\\delta", epsilon: "\\epsilon", zeta: "\\zeta", eta: "\\eta", theta: "\\theta", kappa: "\\kappa", lambda: "\\lambda", mu: "\\mu", nu: "\\nu", xi: "\\xi", rho: "\\rho", sigma: "\\sigma", tau: "\\tau", phi: "\\phi", chi: "\\chi", psi: "\\psi", omega: "\\omega" };
const GREEK_U = { alpha: "α", beta: "β", gamma: "γ", delta: "δ", epsilon: "ε", zeta: "ζ", eta: "η", theta: "θ", kappa: "κ", lambda: "λ", mu: "μ", nu: "ν", xi: "ξ", rho: "ρ", sigma: "σ", tau: "τ", phi: "φ", chi: "χ", psi: "ψ", omega: "ω" };
const WARN_TEX = "#9a4b00";         // --fb-warn, for KaTeX \color
const sub = (n) => "\\mathrm{" + String(n).replace(/_/g, "\\_") + "}";
const sym = (n) => GREEK[n] || sub(n);
// the symbol of a coefficient in text (β), and in TeX in the warn colour when
// the equation does not define it yet
const usym = (n) => GREEK_U[n] || n;
const csym = (n) => coefNames().includes(n) ? sym(n) : `{\\color{${WARN_TEX}}${sym(n)}}`;

// ------------------------------------------------------------ the run's groups

const base = () => D.bases[S.base];
function templateOf(g) {
  if (D.defaults[g]) return g;
  for (const [re, t] of TEMPLATES) if (re.test(g)) return t;
  return null;
}
// Where a group comes from on the current base: the base par file, the
// library, or (a name with a template) new.
function sourceOf(g) {
  if (base().groups[g]) return "base";
  if (D.library.groups[g]) return "library";
  if (templateOf(g)) return "new";
  return null;
}
const sourceDef = (g) => { const s = sourceOf(g); return s === "base" ? base().groups[g] : s === "library" ? D.library.groups[g] : null; };
// a group given by name in the par file or the library (not only as a template)
const named = (g) => ["base", "library"].includes(sourceOf(g));
// The value a key has without the builder: the source group, else the template.
function refOf(g, k) {
  const def = sourceDef(g);
  if (def && def.keys[k]) return { v: def.keys[k].v, src: sourceOf(g) };
  const t = templateOf(g);
  if (t && D.defaults[t].keys[k]) return { v: D.defaults[t].keys[k].v, src: "default" };
  return null;
}
function commentOf(g, k) {
  const def = sourceDef(g), t = templateOf(g);
  const dk = def && def.keys[k], tk = t && D.defaults[t].keys[k];
  return (dk && dk.c) ? dk.c : (tk ? tk.c : "");
}
function schemaKeys(g) {
  const t = templateOf(g), def = sourceDef(g), r = S.g.get(g);
  const order = [];
  if (t) order.push(...D.defaults[t].order);
  if (def) for (const k of def.order) if (!order.includes(k)) order.push(k);
  if (r) for (const k of r.given) if (!order.includes(k)) order.push(k);
  return order;
}
// The record of a group: {k: {key: value}, given: Set(keys written out)},
// created from its source on first use.
function grp(g) {
  if (!S.g.has(g)) {
    if (!sourceOf(g)) return null;
    const k = {}, given = new Set();
    const t = templateOf(g);
    if (t) for (const [key, d] of Object.entries(D.defaults[t].keys)) k[key] = d.v;
    const def = sourceDef(g);
    if (def) for (const [key, d] of Object.entries(def.keys)) { k[key] = d.v; given.add(key); }
    S.g.set(g, { k, given });
  }
  return S.g.get(g);
}
const inRun = (g) => S.g.has(g) || !!sourceOf(g);
function val(g, k, d = undefined) {
  const r = grp(g);
  if (r && k in r.k) return r.k[k];
  const ref = refOf(g, k);
  return ref ? ref.v : d;
}
const asList = (v) => Array.isArray(v) ? v.map(String).filter(s => s !== "") : (v === "" || v == null ? [] : [String(v)]);
const same = (a, b) => JSON.stringify(a) === JSON.stringify(b);
const sameValue = (a, b) => (Array.isArray(a) || Array.isArray(b)) ? same(asList(a), asList(b)) : same(a, b);
function setMany(g, kv) {
  const r = grp(g);
  if (!r) return;
  for (const [k, v] of Object.entries(kv)) { r.k[k] = v; r.given.add(k); }
  render();
}
const set = (g, k, v) => setMany(g, { [k]: v });
function unset(g, k) {
  const r = grp(g);
  if (!r) return;
  const ref = refOf(g, k);
  if (ref) r.k[k] = ref.v; else delete r.k[k];
  r.given.delete(k);
  render();
}
// given, and different from the source (what -p must set)
function changed(g, k) {
  const r = S.g.get(g);
  if (!r || !r.given.has(k)) return false;
  const ref = refOf(g, k);
  return !ref || (ref.src === "default" && sourceOf(g) === "new") || !sameValue(r.k[k], ref.v);
}
// What the builder changed against the sources, as a string: the preset is
// "on" while this is what loading it produced (groups created on first use
// from their source do not count).
const stateKey = () => JSON.stringify([[...S.g].map(([g, r]) => [g, [...r.given].filter(k => changed(g, k)).sort().map(k => [k, r.k[k]])]).filter(x => x[1].length), S.extra]);

// the lists of the equation
const coefNames = () => asList(val("forcing", COEF_LIST_KEY, ""));
const termNames = () => asList(val("forcing", "terms", ""));
const reference = () => String(val("forcing", "reference", "pd"));
const climate = () => String(val("comps", "climate", "snapshots"));
const withSnaps = () => ["snapshots", "rembo"].includes(climate());
// the timeline group ([ctrl] run_step)
const tgroup = () => { const r = String(val("ctrl", "run_step", "ctrl")); return sourceOf(r) ? r : "ctrl"; };
// the named groups of a kind in the par file, the library and the run (prefix stripped)
const namesOf = (prefix) => [...new Set([...Object.keys(base().groups), ...Object.keys(D.library.groups), ...S.g.keys()])].filter(g => g.startsWith(prefix)).map(g => g.slice(prefix.length)).sort();

// ------------------------------------------------------------- import

function loadBase(name) { S.base = name; S.g.clear(); S.extra = []; S.preset = null; S.tsel = null; }
// runme -p tokens into the groups; tokens on groups the builder does not
// know (yelmo physics, regions, ...) are kept and passed through
function parseP(text) {
  const toks = text.match(/(?:[^\s"']+|"[^"]*"|'[^']*')+/g) || [];
  const errs = [];
  for (const tok of toks) {
    if (tok === "-p" || tok === "\\") continue;
    const m = tok.match(/^([\w]+)\.([\w]+)(==|=)(.*)$/s);
    if (!m) { errs.push(`not a key=value: ${tok}`); continue; }
    const [, g, k, op, raw] = m;
    if (!sourceOf(g)) { if (!S.extra.includes(tok)) S.extra.push(tok); continue; }
    const ref = refOf(g, k);
    const items = (op === "==" ? raw.split(",") : [raw]).map(s => s.replace(/^["']|["']$/g, ""));
    const refel = ref ? (Array.isArray(ref.v) ? (ref.v[0] ?? "") : ref.v) : "";
    const conv = (s) => typeof refel === "boolean" ? /^(true|t|1|\.true\.)$/i.test(s) : typeof refel === "number" ? num(s) : s;
    const v = (op === "==" || (ref && Array.isArray(ref.v))) ? items.map(conv) : conv(items[0]);
    const r = grp(g);
    r.k[k] = v; r.given.add(k);
  }
  return errs;
}
function loadPreset(p) {
  loadBase(p.base);
  parseP(p.p.join(" "));
  S.preset = p.name; S.presetKey = stateKey();
  render();
}
const presetState = () => S.preset && S.presetKey === stateKey() ? "on" : "mod";

// ------------------------------------------------------------------- helpers

function h(tag, attrs = {}, ...children) {
  const el = document.createElement(tag);
  for (const [a, v] of Object.entries(attrs)) {
    if (a === "class") el.className = v;
    else if (a === "html") el.innerHTML = v;
    else if (a.startsWith("on")) el.addEventListener(a.slice(2), v);
    else if (v !== null && v !== undefined && v !== false) el.setAttribute(a, v === true ? "" : v);
  }
  for (const c of children.flat()) if (c != null) el.append(c.nodeType ? c : document.createTextNode(String(c)));
  return el;
}
const fmtNum = (x) => {
  if (typeof x !== "number") return String(x);
  if (Number.isInteger(x) && Math.abs(x) >= 1000 && x % 1000 === 0) return (x / 1000) + "e3";
  return Number.isInteger(x) ? String(x) : String(+x.toPrecision(6));
};
const fmtK = (x) => fmtNum(x).replace(/e3$/, "\\,\\mathrm{k}");
const fmtV = (x) => Number.isFinite(x) ? fmtNum(+x.toPrecision(4)) : "–";
function tex(s, display = true) {
  return h("div", { class: display ? "fb-eq" : "fb-eq-inline", html: katex.renderToString(s, { displayMode: display, throwOnError: false }) });
}
const texi = (s) => h("span", { class: "fb-sym", html: katex.renderToString(s, { throwOnError: false }) });
const num = (s) => { const x = Number(String(s).replace(/d/i, "e")); return Number.isFinite(x) ? x : NaN; };

// Known option lists, by "group.key" (or "template.key"); otherwise from the
// first line of the key's comment (a|b|c).
const OPTIONS = {
  "ctrl.tstep_method": ["const", "rel", "cal"], "spinup.tstep_method": ["const", "rel", "cal"], "transient.tstep_method": ["const", "rel", "cal"],
  "barysealevel.method": ["const", "file", "fastiso"],
  "orbit.method": ["forcing_time", "fixed_time", "fixed_orbit"], "ice_bnd.method": ["file", "scaled"],
  [COEF_T + ".from"]: ["", "air_anomaly"],
  "term.pattern": ["snapshots", "uniform"],
  "taper.outside": ["ref", "mean"],
};

// The help under a field, in plain words, by "group.key" or "template.key";
// consulted before the namelist comment. The time keys are written for ctrl
// and serve the timeline group of the run ([ctrl] run_step).
const HELP = {
  "ctrl.tstep_method": "How the forcing time τ relates to the model time: held at tstep_const (const), equal to the model time (rel), or calendar years (cal).",
  "ctrl.tstep_const": "The forcing time τ [yr] the whole run is held at (const); past < 0, e.g. -21000 for the LGM.",
  "ctrl.time_init": "The model time [yr] at the start of the run.",
  "ctrl.time_end": "The model time [yr] at the end of the run.",
  "ctrl.calendar_ref": "The calendar year [yr CE] where τ = 0 (cal): τ = t − calendar_ref.",

  "forcing.reference": "The reference snapshot S₀: the climate every term adds to.",
  "snapshots.lapse": "The lapse rates Γ [K/m], annual and summer: every snapshot is reduced to sea level with them, and the climate brought back to the model surface.",
  "snapshots.f_p": "Precipitation sensitivity: P scales by exp(β_p ΔT) with elevation and in uniform terms.",

  ["term." + TERM_COEF_KEY]: "The coefficient that scales this term; empty = the constant 1 (the full anomaly of a reconstruction).",
  "term.pattern": "snapshots: coefficient × (snap1 − snap0); uniform: coefficient × a gain, the same everywhere.",
  "term.snap0": "The snapshot where the coefficient is 0 (e.g. piControl). Empty: snap1 is itself an anomaly or ratio.",
  "term.snap1": "The snapshot where the coefficient is 1 (e.g. LGM).",
  "term.fields": "The fields the term changes. tas, to, so add; pr multiplies (fixed rule).",
  "term.f_ta": "Gain [K] on the air temperature: ΔT = f_ta × the coefficient, everywhere; pr follows it through β_p.",
  "term.f_to": "Gain [K] on the ocean temperature: ΔΘ = f_to × the coefficient, at every depth.",
  "term.f_so": "Gain [psu] on the salinity: ΔS = f_so × the coefficient, at every depth.",

  [COEF_T + ".method"]: "How the coefficient moves in time: the curve of a ramp or wave, or the rule of a controller.",
  [COEF_T + ".units"]: "The units of the coefficient, a label in the 1D output.",
  [COEF_T + ".f_min"]: "The lower limit of the coefficient; a rising ramp starts here.",
  [COEF_T + ".f_max"]: "The upper limit of the coefficient; a falling ramp starts here.",
  [COEF_T + ".df_dt_max"]: "The rate [1/yr] of the ramp; for a controller, the fastest it may move.",
  [COEF_T + ".df_sign"]: "+1: from f_min upwards; -1: from f_max downwards.",
  [COEF_T + ".dt_init"]: "Years [yr] at the start where the coefficient is held at its start value.",
  [COEF_T + ".dt_ramp"]: "The length [yr] of the ramp, or the period of the wave.",
  [COEF_T + ".dt_conv"]: "The years [yr] after the ramp to move to f_conv.",
  [COEF_T + ".f_conv"]: "The value held after the ramp and the convergence.",
  [COEF_T + ".dt_ave"]: "The window [yr] over which the controller averages the ice-volume change.",
  [COEF_T + ".tol"]: "The ice-volume change rate the controller treats as equilibrium.",
  [COEF_T + ".with_kill"]: "Stop the run once the ice sheet is in equilibrium at a limit (f_min or f_max).",
  [COEF_T + ".from"]: "What the model derives the coefficient from: the mean air temperature anomaly (air_anomaly).",
  [COEF_T + ".factor"]: "The coefficient = factor × the mean sea-level air temperature anomaly.",
  [COEF_T + ".region"]: "Where the anomaly is averaged: all = the whole climate grid, or a region, e.g. zone:open_ocean.",
  [COEF_T + ".anomaly"]: "On: the box temperature minus its value at the start; off: the temperature itself.",
  [COEF_T + ".series_file"]: "The file of the record. Not in the list? Any two-column table or netCDF time series in input/ works for the model; the list is what scripts/forcing-builder/export.jl found when the page was built. Put the file in input/ and re-run the export.",
  [COEF_T + ".series_var"]: "The variable of the netCDF file that holds the record.",
  [COEF_T + ".series_time_name"]: "The name of the time variable in the netCDF file.",
  [COEF_T + ".norm_t0"]: "Time where the normalised record is 0; two times = the mean over that window.",
  [COEF_T + ".norm_t1"]: "Time where the normalised record is 1; two times = the mean over that window.",

  "snapshot.tas": "The clim_db group(s) of the air temperature: one monthly group, or two for annual + summer temperature.",
  "snapshot.pr": "The clim_db group(s) of the precipitation: one, or two summed (snowfall + rainfall).",
  "snapshot.zs": "The clim_db group of the surface elevation of the climate data, for the reduction to sea level.",
  "snapshot.to": "The clim_db group of the ocean temperature profiles.",
  "snapshot.so": "The clim_db group of the ocean salinity profiles.",
  "snapshot.exp": "The experiment name filled into {exp} of the file names (lgm, piControl, ...).",
  "snapshot.monthly": "On: tas and pr are monthly fields; off: tas is annual + summer, pr annual.",
  "snapshot.clim_sigma": "Smoothing [km] of the atmosphere fields when they are read; 0 = none.",
  "snapshot.ocn_sigma": "Smoothing [km] of the ocean fields when they are read; 0 = none.",

  "esm.experiment": "The ESM experiment of the projection (ssp585, 1pctCO2, ...); ctrl = the reference climate only.",
  "esm.esm_name": "The ESM whose anomalies are read.",
  "esm.use_esm": "Read the ESM anomalies; off: the reference climate, or the homogeneous warming of 1pctCO2.",
  "esm.use_hist": "Read the historical period of the ESM.",
  "esm.use_proj": "Read the projection period of the ESM.",

  "taper.active": "Pull the climate towards the reference where the weight is below 1.",
  "taper.weight": "The clim_db group of the weight w: 1 keeps the composed climate, 0 the reference.",
  "taper.outside": "Where w < 1, towards the reference (ref), or the reference plus the mean anomaly where w = 1 (mean).",
  "bias.active": "Add a temperature correction and multiply the precipitation by a correction factor.",
  "bias.dtas": "The clim_db group of the temperature correction [K].",
  "bias.fpr": "The clim_db group of the precipitation factor [1].",
  "bias.f_tas": "The weight of the temperature correction.",
  "bias.f_pr": "The weight of the precipitation factor (applied to factor − 1).",
  "orbit.method": "The orbit of the insolation: at the forcing time τ (forcing_time), at one fixed time (fixed_time), or given elements (fixed_orbit).",
  "orbit.time": "The time [yr] of the fixed orbit, relative to 1950 (past < 0).",
  "orbit.ecc": "The eccentricity of the fixed orbit.",
  "orbit.obl": "The obliquity [deg] of the fixed orbit.",
  "orbit.per": "The longitude of perihelion [deg] of the fixed orbit.",
  "barysealevel.method": "The sea level: a constant (const), a record read at τ (file), or from the ice volume of the FastIsostasy domains (fastiso).",
  "barysealevel.bsl_init": "The constant sea level [m] relative to present day.",
  "barysealevel.sl_path": "The file of the sea-level record, read at τ.",
  "barysealevel.sl_name": "The variable of a netCDF sea-level file; none for a table.",
  "ice_bnd.active": "Prescribe the ice thickness in the fixed and relaxation cells at the domain border.",
  "ice_bnd.method": "file: the ice thickness at τ from the file; scaled: a coefficient × one slice of it.",
  "ice_bnd.H_ice": "The clim_db group of the prescribed ice thickness [m].",
  "ice_bnd.coef": "scaled: the coefficient that scales the slice.",
  "ice_bnd.slice_time": "scaled: the time [yr] of the slice.",
};
// the help of a key: HELP, else the namelist comment in full
function helpOf(g, k) {
  const t = templateOf(g);
  for (const hk of [g + "." + k, t && t + "." + k, g === tgroup() ? "ctrl." + k : null]) if (hk && HELP[hk]) return HELP[hk];
  return (commentOf(g, k) || "").split("\n").map(s => s.trim()).join(" ");
}
function optionsOf(g, k) {
  const t = templateOf(g);
  if (OPTIONS[g + "." + k]) return OPTIONS[g + "." + k];
  if (t && OPTIONS[t + "." + k]) return OPTIONS[t + "." + k];
  const m = (commentOf(g, k) || "").split("\n")[0].match(/^([\w.-]+(\|[\w.-]+)+)/);
  if (m) return m[1].split("|");
  return null;
}

// A labelled widget of one key: checkbox, select, number or text; a badge
// when the value differs from the source.
// opts: {label (text or a node), options, labels: {value: label}, help,
// chips: the clim_db groups to pick from (a chip list), single}
function field(g, k, opts = {}) {
  const r = grp(g);
  if (!r) return h("div", { class: "fb-field fb-missing" }, `[${g}]: no such group`);
  const v = val(g, k);
  if (v === undefined) return h("div", { class: "fb-field fb-missing" }, `${g}.${k}: not in the schema`);
  if (opts.chips) return chipList(g, k, opts.chips, opts);
  let input;
  const options = opts.options || optionsOf(g, k);
  const labels = opts.labels || {};
  if (typeof v === "boolean") {
    input = h("input", { type: "checkbox", onchange: (e) => set(g, k, e.target.checked) });
    input.checked = v;
  } else if (options && !Array.isArray(v)) {
    input = h("select", { onchange: (e) => set(g, k, e.target.value) },
      ...options.map(o => h("option", { value: o }, labels[o] ?? o)), ...(options.includes(String(v)) ? [] : [h("option", { value: v }, String(v))]));
    input.value = String(v);
  } else if (typeof v === "number") {
    input = h("input", { type: "text", class: "fb-num", value: fmtNum(v), onchange: (e) => { const x = num(e.target.value); if (Number.isFinite(x)) set(g, k, x); else render(); } });
  } else if (Array.isArray(v)) {
    const isnum = v.every(x => typeof x === "number");
    input = h("input", { type: "text", class: "fb-text", value: v.map(fmtNum).join(", "),
      onchange: (e) => { const items = e.target.value.split(/[,\s]+/).filter(s => s !== ""); set(g, k, isnum ? items.map(num) : items); } });
  } else {
    input = h("input", { type: "text", class: "fb-text", value: String(v), onchange: (e) => set(g, k, e.target.value) });
  }
  return fieldRow(g, k, input, opts);
}
// One row: the label, the widget with its badge, and the help in full on a
// second line under the widget. The
// badge: "new" for a key of a new group, "changed" where -p must set it.
function fieldRow(g, k, input, opts = {}) {
  const ch = changed(g, k), ref = refOf(g, k), src = sourceOf(g);
  const was = ref && ref.src !== "default" ? ` (there: ${fmtNml(ref.v)})` : "";
  const badge = !ch ? h("span", {}) : src === "new"
    ? h("span", { class: "fb-src fb-src-new", title: "a key of a new group: written in the block to add" }, "new")
    : h("span", { class: "fb-src fb-src-override", title: `differs from the ${src === "library" ? "library" : "par file"}${was}; set with -p` }, "changed");
  const comment = commentOf(g, k);
  const help = opts.help !== undefined ? opts.help : helpOf(g, k);
  const label = opts.label && opts.label.nodeType ? opts.label : h("code", {}, opts.label || k);
  const tip = comment || (help ? "" : helpOf(g, k));
  return h("div", { class: "fb-field" },
    h("label", { title: tip ? `${g}.${k}: ${tip}` : `${g}.${k}` }, label), h("span", { class: "fb-val" }, input, badge), help ? h("span", { class: "fb-help" }, help) : null);
}
function fields(g, keys, opts = {}) { return h("div", { class: "fb-fields" }, ...keys.map(k => field(g, k, opts[k] || {}))); }
// A vector-of-strings key as checkboxes over a fixed set (the fields of a term).
function checkList(g, k, options, marks = {}) {
  const cur = asList(val(g, k, ""));
  const input = h("span", { class: "fb-checks" }, ...options.map(o => {
    const cb = h("input", { type: "checkbox", onchange: (e) => set(g, k, options.filter(x => x === o ? e.target.checked : cur.includes(x))) });
    cb.checked = cur.includes(o);
    return h("label", { title: marks[o] || "" }, cb, o + (marks[o] ? " ✗" : ""));
  }));
  return fieldRow(g, k, input);
}
// A name key as a select over known names, with "new name…": choosing it
// shows a text input next to the select; Enter or leaving it sets the key.
// items: [[value, label]]; the current value is added when not among them.
const NEW_NAME = "\u0000new";
function namePicker(g, k, items, opts = {}) {
  const v = String(val(g, k, ""));
  const all = items.some(([x]) => x === v) ? items : [...items, [v, v === "" ? "(none)" : `${v} (not defined)`]];
  const sel = h("select", {}, ...all.map(([x, lab]) => h("option", { value: x }, lab)), h("option", { value: NEW_NAME }, "new name…"));
  sel.value = v;
  const txt = h("input", { type: "text", class: "fb-num fb-newname", placeholder: "a name", hidden: true });
  let done = false;
  const commit = () => {
    if (done) return;
    done = true;
    const n = txt.value.trim();
    if (/^\w+$/.test(n)) set(g, k, n); else render();
  };
  sel.addEventListener("change", () => {
    if (sel.value === NEW_NAME) { txt.hidden = false; txt.focus(); } else set(g, k, sel.value);
  });
  txt.addEventListener("keydown", (e) => { if (e.key === "Enter") commit(); else if (e.key === "Escape") { done = true; render(); } });
  txt.addEventListener("blur", commit);
  return fieldRow(g, k, h("span", { class: "fb-picker" }, sel, txt), opts);
}
// The snapshot of a reference or a term: the snapshots of the par file, the
// library and those the run created new; snap0 may be none.
function snapPicker(g, k, opts = {}) {
  const P = "snapshot_";
  const where = (s) => ({ base: "in par file", library: "in library", new: "new" })[sourceOf(P + s)] || "?";
  // the current name too: a new one is created by the snapshot block, after this picker
  const v = String(val(g, k, "")), names = namesOf(P);
  if (v && !names.includes(v) && sourceOf(P + v)) names.push(v);
  const items = names.map(s => [s, `${s} (${where(s)})`]);
  if (opts.none) items.unshift(["", opts.none]);
  return namePicker(g, k, items, opts);
}
// The coefficient of a term: the defined ones, the constant 1, or a new name.
function coefPicker(g, k, opts = {}) {
  return namePicker(g, k, [["", "1: no coefficient"], ...coefNames().map(n => [n, n])], opts);
}
// The clim_db groups of a key as chips (× removes one) and a select that
// adds one; single: one group, adding replaces it.
function chipList(g, k, options, opts = {}) {
  const v = val(g, k, "");
  const cur = asList(v);
  const put = (list) => set(g, k, opts.single ? (list[0] ?? "") : list);
  const chips = cur.map((c, i) => h("span", { class: "fb-chip" + (options.includes(c) ? "" : " fb-chip-bad"), title: options.includes(c) ? c : `${c}: not a group of input/clim_db.nml` }, c,
    h("button", { class: "fb-chip-x", title: `remove ${c}`, onclick: () => put(cur.filter((_, j) => j !== i)) }, "×")));
  const add = h("select", { class: "fb-chip-add", onchange: (e) => { if (e.target.value) put(opts.single ? [e.target.value] : [...cur, e.target.value]); } },
    h("option", { value: "" }, opts.single && cur.length ? "replace with…" : "add a group…"), ...options.map(o => h("option", { value: o }, o)));
  return fieldRow(g, k, h("span", { class: "fb-chips" }, ...chips, add), opts);
}
const climdbNames = () => Object.keys(D.climdb).sort();
// The keys a group sets that the cards do not show (the rare ones), as a line.
function alsoSets(g, shown) {
  const r = S.g.get(g);
  if (!r) return null;
  const tk = (D.defaults[templateOf(g)] || { keys: {} }).keys;
  const rest = [...r.given].filter(k => !shown.includes(k) && (!(k in tk) || !sameValue(r.k[k], tk[k].v)));
  if (!rest.length) return null;
  return h("div", { class: "fb-vals" }, "also set: " + rest.map(k => `${k} = ${fmtNml(r.k[k])}`).join(", ") + " (in the par file; not edited here)");
}
// The gist of a group's comment for a list: its first clause (up to the first
// " (", ": " or ". " past a few words), at most 70 characters.
function gist(c) {
  const s = String(c || "").split("\n").map(x => x.trim()).join(" ").trim();
  if (!s) return "";
  // not at "et al." nor at a year in parentheses (a citation)
  const cut = [...s.matchAll(/ \((?!\d)|: |(?<!\bal)\. /g)].map(m => m.index).filter(i => i >= 15);
  const t = (cut.length ? s.slice(0, Math.min(...cut)) : s).replace(/[.:]$/, "");
  return t.length > 70 ? t.slice(0, 69).trimEnd() + "…" : t;
}
// A key of a library group as a string ("" when the library does not set it).
const libKey = (g, k) => { const d = D.library.groups[g]; return d && d.keys[k] ? String(d.keys[k].v) : ""; };
// The library groups of the yelmox_bipolar driver (the ocean box model): not
// offered here, documented on its page.
const bipolarOnly = (g) => libKey(g, "from").startsWith("obm_") || libKey(g, "pattern") === "replace";
const groupGist = (g) => { const s = sourceOf(g); return gist(s === "base" ? base().groups[g].c : s === "library" ? D.library.groups[g].c : ""); };

// The add box under a list of names (the coefficients, the terms of
// [forcing]): a named group of the par file or the library, each with the
// gist of its comment, or a new name (from the template); not the groups of
// the yelmox_bipolar driver.
function addBox(k, prefix, what, onnew) {
  const list = asList(val("forcing", k, ""));
  const cands = namesOf(prefix).filter(c => !list.includes(c) && !bipolarOnly(prefix + c));
  const where = (c) => sourceOf(prefix + c) === "library" ? " (in library)" : sourceOf(prefix + c) === "base" ? " (in par file)" : "";
  const sel = h("select", {}, h("option", { value: "" }, "one of the par file or the library…"),
    ...cands.map(c => { const gs = groupGist(prefix + c); return h("option", { value: c, title: prefix + c }, c + (gs ? " — " + gs : "") + where(c)); }));
  const txt = h("input", { type: "text", class: "fb-num", placeholder: "or a new name" });
  const add = () => {
    const n = sel.value || txt.value.trim();
    if (!n || list.includes(n) || !/^\w+$/.test(n)) return;
    const isnew = !named(prefix + n) && !S.g.has(prefix + n);
    grp("forcing"); grp(prefix + n);
    if (isnew && onnew) onnew(prefix + n);
    set("forcing", k, [...list, n]);
  };
  return h("div", { class: "fb-card fb-card-add" }, h("div", { class: "fb-names" }, h("span", {}, `add a ${what}:`), sel, txt, h("button", { class: "fb-btn", onclick: add }, "add")));
}
function removeBtn(k, name) {
  return h("button", { class: "fb-btn fb-remove", title: `remove ${name} from [forcing] ${k}`, onclick: () => set("forcing", k, asList(val("forcing", k, "")).filter(x => x !== name)) }, "remove");
}

// ------------------------------------------------------------------ timeline

function timeline() {
  const T = tgroup();
  const method = String(val(T, "tstep_method", "const"));
  const t0 = num(val(T, "time_init", 0)), t1 = num(val(T, "time_end", 0));
  const tc = num(val(T, "tstep_const", 0));
  const cref = num(val("ctrl", "calendar_ref", 1950));
  const n = 400;
  const t = Array.from({ length: n + 1 }, (_, i) => t0 + (t1 - t0) * i / n);
  const tf = t.map(x => method === "const" ? tc : (method === "cal" ? x - cref : x));
  return { T, method, t0, t1, tc, cref, t, tf, tf0: tf[0], tfmin: Math.min(...tf), tfmax: Math.max(...tf) };
}
function evalTime(tl) {
  const lo = Math.min(tl.t0, tl.t1), hi = Math.max(tl.t0, tl.t1);
  return S.tsel === null || !Number.isFinite(S.tsel) ? tl.t0 : Math.min(Math.max(S.tsel, lo), hi);
}
const tauAt = (tl, t) => tl.method === "const" ? tl.tc : (tl.method === "cal" ? t - tl.cref : t);
const atT = (tl, f, t) => f ? interp({ time: tl.t, f }, t) : NaN;

// ----------------------------------------------------------------- tsgen (JS)

const KINDS = { record: ["series"], constant: ["const"], ramp: ["ramp-slope", "ramp-time", "ramp-time-step", "sin"], controller: ["exp", "PI42", "H312b", "H312PID", "H321PID", "PID1"], derived: ["derived"] };
const KIND_LABEL = { record: "series (a time series from a file)", constant: "constant", ramp: "ramp or wave", controller: "controller (feedback)", derived: "derived (set by the model)" };
// the kind of a coefficient: tsgen "const" with a rate is a ramp
function kindOf(p) {
  const m = String(p.method);
  if (m === "const" && num(p.df_dt_max) > 0) return "ramp";
  return Object.keys(KINDS).find(k => KINDS[k].includes(m)) || "constant";
}
function seriesOf(file, varname) {
  const s = D.series[file];
  if (!s) return null;
  if (s.kind === "netcdf") { const v = s.vars[varname]; return v ? { time: v.time, f: v.f, header: `${file}:${varname}`, n: v.n } : null; }
  return s;
}
function interp(ser, t) {
  const T = ser.time, F = ser.f, n = T.length;
  if (!n) return NaN;
  if (t <= T[0]) return F[0]; if (t >= T[n - 1]) return F[n - 1];
  let lo = 0, hi = n - 1;
  while (hi - lo > 1) { const m = (lo + hi) >> 1; if (T[m] <= t) lo = m; else hi = m; }
  const w = (t - T[lo]) / (T[hi] - T[lo]);
  return F[lo] * (1 - w) + F[hi] * w;
}
function seriesMean(ser, tw) {
  if (tw.length === 1) return interp(ser, tw[0]);
  const lo = Math.min(...tw), hi = Math.max(...tw);
  const sel = ser.f.filter((f, i) => ser.time[i] >= lo && ser.time[i] <= hi);
  return sel.length ? sel.reduce((a, b) => a + b, 0) / sel.length : interp(ser, (lo + hi) / 2);
}
// The parameters of [coef_<name>]; normalisation is active when the group
// itself sets norm_t0 and norm_t1 (the model's rule).
function coefPar(name) {
  const g = COEF_PREFIX + name, r = grp(g);
  if (!r) return null;
  const p = {};
  for (const k of schemaKeys(g)) p[k] = val(g, k);
  p.norm = r.given.has("norm_t0") && r.given.has("norm_t1");
  return p;
}
function seriesFile(p) {
  if (p.method !== "series" || !p.series_file) return null;
  let ser = seriesOf(p.series_file, p.series_var);
  if (!ser) return null;
  if (p.norm) {
    const v0 = seriesMean(ser, asList(p.norm_t0).map(num)), v1 = seriesMean(ser, asList(p.norm_t1).map(num));
    ser = { ...ser, f: ser.f.map(f => (f - v0) / (v1 - v0)), norm: [v0, v1] };
  }
  return ser;
}
function evalOpen(p, te) {
  const fmin = num(p.f_min), fmax = num(p.f_max), s = num(p.df_sign) >= 0 ? 1 : -1;
  const tau = Math.max(0, te - num(p.dt_init));
  const fstart = s > 0 ? fmin : fmax;
  const clamp = (x) => Math.min(Math.max(x, fmin), fmax);
  switch (p.method) {
    case "const": case "ramp-slope": return clamp(fstart + s * num(p.df_dt_max) * tau);
    case "ramp-time": return clamp(fstart + s * Math.abs(fmax - fmin) / num(p.dt_ramp) * tau);
    case "ramp-time-step": {
      const dr = num(p.dt_ramp), dc = num(p.dt_conv), fc = num(p.f_conv), fext = s > 0 ? fmax : fmin;
      if (tau <= dr) return fstart + (fext - fstart) * (tau / dr);
      if (tau <= dr + dc) return fext + (fc - fext) * ((tau - dr) / dc);
      return fc;
    }
    case "sin": { const T = te <= num(p.dt_init) ? 0 : te; return 0.5 * (fmax + fmin) + 0.5 * (fmax - fmin) * Math.sin(2 * Math.PI * T / num(p.dt_ramp)); }
    default: return fstart;
  }
}
function applyWindow(p, t, f) {
  const t0 = num(p.win_t0), t1 = num(p.win_t1);
  if (t1 <= t0 || t <= t0 || t >= t1) return f;
  const ref = num(p.win_ref);
  if (p.win_above && f <= ref) return f;
  const taper = num(p.win_taper);
  const w = taper > 0 ? Math.min(1, (t - t0) / taper, (t1 - t) / taper) : 1;
  return f + w * (num(p.win_fac) - 1) * (f - ref);
}
// The coefficient on the run's timeline: {p, ser, f, kind}; derived: f = null.
function coefOnRun(name, tl) {
  const p = coefPar(name);
  if (!p) return { p: null, ser: null, f: null, kind: null, missing: true };
  const kind = kindOf(p);
  if (kind === "derived") return { p, ser: null, f: null, kind, missing: false };
  const ser = seriesFile(p);
  const f = tl.tf.map(tf => applyWindow(p, tf, p.method === "series" ? (ser ? interp(ser, tf) : NaN) : evalOpen(p, tf - tl.tf0)));
  return { p, ser, f, kind, missing: false };
}
function termWeight(t, tl) {
  if (t.coef === "") return tl.t.map(() => 1);
  if (!coefNames().includes(t.coef)) return null;
  const s = coefOnRun(t.coef, tl);
  return s.f && s.f.some(Number.isFinite) ? s.f : null;
}

// ------------------------------------------------------------------- figures

const WR = 700, HR = 200, WL = 560, HL = 120;
function linePlot(series, { x = "t [yr]", y = "", rules = [], title = "", vrules = [], band = null, width = WR, height = HR, tmark = null }) {
  const data = series.flatMap(s => s.t.map((t, i) => ({ t, f: s.f[i], label: s.label })));
  const marks = [Plot.ruleY([0], { stroke: "#999" })];
  if (band) marks.push(Plot.rectX([band], { x1: "x0", x2: "x1", fill: "#ffd", stroke: "none" }));
  marks.push(Plot.line(data, { x: "t", y: "f", stroke: "label", strokeWidth: 1.5 }));
  for (const r of rules) marks.push(Plot.ruleY([r.y], { stroke: r.stroke || "#c33", strokeDasharray: "4 3" }), Plot.text([{ t: data.length ? data[0].t : 0, f: r.y }], { x: "t", y: "f", text: () => r.label, textAnchor: "start", dy: -6, fill: r.stroke || "#c33" }));
  for (const r of vrules) marks.push(Plot.ruleX([r.x], { stroke: r.stroke || "#36c", strokeDasharray: "4 3" }), Plot.text([{ t: r.x }], { x: "t", text: () => r.label, frameAnchor: "top", dy: 4, fill: r.stroke || "#36c" }));
  if (tmark !== null) marks.push(Plot.ruleX([tmark], { stroke: "#555", strokeWidth: 1 }));
  const ys = data.map(d => d.f).filter(Number.isFinite).concat(rules.map(r => r.y));
  let domain;
  if (ys.length) { const lo = Math.min(...ys), hi = Math.max(...ys); if (hi - lo < 1e-9) domain = [lo - 1, hi + 1]; }
  const el = Plot.plot({ width, height, marginLeft: 55, title: title || undefined, x: { label: x }, y: { label: y, grid: true, domain },
    color: series.length > 1 ? { legend: true } : undefined, marks });
  return h("figure", { class: "fb-fig" }, el);
}
const smallPlot = (series, opts = {}) => linePlot(series, { width: WL, height: HL, ...opts });

// ------------------------------------------------------------------- layout

function section(id, title, hint, ...content) {
  return h("section", { class: "fb-section", id: "fb-" + id }, h("h2", {}, title, hint ? h("span", { class: "fb-hint" }, hint) : null), ...content);
}
function eqblock(title, ...content) { return h("div", { class: "fb-eqblock" }, h("h3", {}, title), ...content); }
function details(id, title, ...content) {
  const d = h("details", { class: "fb-details" }, h("summary", {}, title), ...content);
  if (S.open.has(id)) d.open = true;
  d.addEventListener("toggle", () => { if (d.open) S.open.add(id); else S.open.delete(id); });
  return d;
}
function note(text, cls = "fb-note") { return h("p", { class: cls }, text); }
function warn(text) { return note(text, "fb-warn"); }
function vals(text) { return h("div", { class: "fb-vals" }, text); }
// The hints where a piece is missing (a clim_db group, a group in no file yet),
// and the note on what belongs to the yelmox_bipolar driver only.
function climdbHint() {
  return h("p", { class: "fb-note" }, "Not in the list? Add a group to input/clim_db.nml (file, variable, units, grid, time axis; the format: ",
    h("a", { href: "database-namelists.html" }, "database namelists"), "), then re-run the export. The file itself is a FesmData product on the climate grid.");
}
const newHint = (g) => sourceOf(g) !== "new" ? null
  : note("A new group is not yet in any file: the written-out run (right column, par-file tab) gives its block; add it to the par file for this setup or to input/forcing_db.nml to share it.");
function bipolarNote(text) {
  return h("p", { class: "fb-note" }, text + "; used by the yelmox_bipolar driver only, see ", h("a", { href: "yelmox-bipolar.html" }, "yelmox_bipolar"), ".");
}
// Where a group is defined: the badge and its tooltip.
const SRC_BADGE = {
  base: ["in par file", () => `defined in ${base().file}`],
  library: ["in library", () => `defined in ${D.library.files.join(", ")}, shared by every par file`],
  new: ["new", () => "not defined in the par file or the library yet: the written-out run gives the block to add"],
};
const srcBadge = (g) => { const s = sourceOf(g), b = SRC_BADGE[s]; return b ? h("span", { class: "fb-src fb-src-" + s, title: b[1]() }, b[0]) : null; };
// The legend of the badges, once under the equation.
function badgeLegend() {
  const item = (cls, text, what) => h("span", {}, h("span", { class: "fb-src " + cls }, text), " " + what);
  return h("p", { class: "fb-legend" }, item("fb-src-base", "in par file", base().file.replace(/^yelmox\//, "")), item("fb-src-library", "in library", "forcing_db.nml, shared"),
    item("fb-src-new", "new", "not in any file yet; the written-out run gives the block to add"), item("fb-src-override", "changed", "set with -p"));
}
function card(head, g, ...content) {
  return h("div", { class: "fb-card" }, h("div", { class: "fb-cardhead" }, h("code", {}, head), g ? srcBadge(g) : null), ...content);
}

// ------------------------------------------------------------------ top bar

function topBar() {
  const state = presetState();
  const presets = (D.presets || []).map(p => h("button", { class: "fb-preset" + (S.preset === p.name ? (state === "on" ? " on" : " mod") : ""), title: `${D.bases[p.base].file} -p ${p.p.join(" ")}`, onclick: () => loadPreset(p) }, p.name));
  const sel = h("select", {}, ...Object.keys(D.bases).map(b => h("option", { value: b }, b)));
  sel.value = S.base;
  return h("div", { class: "fb-top" },
    h("div", { class: "fb-toprow" }, h("span", { class: "fb-hint" }, "start from a case of the docs figures:"), h("div", { class: "fb-presets" }, ...presets),
      S.preset ? h("span", { class: "fb-help" }, state === "on" ? S.preset : `${S.preset}, modified`) : null),
    h("div", { class: "fb-toprow" }, h("span", { class: "fb-hint" }, "or from a base par file as it is:"), sel,
      h("button", { class: "fb-btn", onclick: () => { loadBase(sel.value); render(); } }, "load"),
      h("span", { class: "fb-help" }, `now: ${base().file}, domain ${val("domain", "name", "?")}, climate backend ${climate()}`)));
}

// ------------------------------------------------------------- left: time

function secTime(tl) {
  const T = tl.T;
  const keys = ["tstep_method", ...(tl.method === "const" ? ["tstep_const"] : []), "time_init", "time_end"];
  return section("time", "Forcing time", `τ, set by [${T}] tstep_method`,
    h("div", { class: "fb-card" },
      fields(T, keys, { tstep_method: { labels: { const: "const: τ = tstep_const", rel: "rel: τ = t", cal: "cal: τ = t − calendar_ref" } } }),
      tl.method === "cal" ? fields("ctrl", ["calendar_ref"]) : null,
      note(tl.method === "const" ? `The forcing is held at τ = ${fmtNum(tl.tc)} yr while the model time runs from ${fmtNum(tl.t0)} to ${fmtNum(tl.t1)} yr: every coefficient, time-dependent snapshot and record is read once, at τ.`
        : tl.method === "rel" ? "The forcing follows the model time: the coefficients and the records are read at τ = t (paleo runs)."
        : `Calendar years: τ = t − ${fmtNum(tl.cref)} (ESM runs).`)));
}

// --------------------------------------------------------- left: equation

const METHOD_EQ = {
  "const": (p, a) => `${a} = \\min\\!\\big(\\max(f_{\\min} + s\\,\\dot f_{\\max}\\,\\tau_e,\\ f_{\\min}),\\ f_{\\max}\\big)`,
  "ramp-slope": (p, a) => `${a} = \\min\\!\\big(\\max(f_{\\mathrm{start}} + s\\,\\dot f_{\\max}\\,\\tau_e,\\ f_{\\min}),\\ f_{\\max}\\big)`,
  "ramp-time": (p, a) => `${a} = \\min\\!\\big(\\max(f_{\\mathrm{start}} + s\\,\\tfrac{|f_{\\max}-f_{\\min}|}{\\Delta t_{\\mathrm{ramp}}}\\,\\tau_e,\\ f_{\\min}),\\ f_{\\max}\\big)`,
  "ramp-time-step": (p, a) => `${a} = \\begin{cases} f_{\\mathrm{start}} + (f_{\\mathrm{ext}}-f_{\\mathrm{start}})\\,\\tau_e/\\Delta t_{\\mathrm{ramp}} & \\tau_e \\le \\Delta t_{\\mathrm{ramp}} \\\\ f_{\\mathrm{ext}} + (f_{\\mathrm{conv}}-f_{\\mathrm{ext}})\\,(\\tau_e-\\Delta t_{\\mathrm{ramp}})/\\Delta t_{\\mathrm{conv}} & \\le \\Delta t_{\\mathrm{ramp}}+\\Delta t_{\\mathrm{conv}} \\\\ f_{\\mathrm{conv}} & \\text{after} \\end{cases}`,
  "sin": (p, a) => `${a} = \\tfrac{f_{\\max}+f_{\\min}}{2} + \\tfrac{f_{\\max}-f_{\\min}}{2}\\,\\sin\\!\\big(2\\pi\\, t_e/\\Delta t_{\\mathrm{ramp}}\\big)`,
};
const METHOD_KEYS = {
  "const": ["f_min", "f_max", "df_dt_max", "df_sign", "dt_init"],
  "ramp-slope": ["f_min", "f_max", "df_dt_max", "df_sign", "dt_init"],
  "ramp-time": ["f_min", "f_max", "dt_ramp", "df_sign", "dt_init"],
  "ramp-time-step": ["f_min", "f_max", "dt_ramp", "dt_conv", "f_conv", "df_sign", "dt_init"],
  "sin": ["f_min", "f_max", "dt_ramp", "dt_init"],
};
const CONTROLLER_KEYS = ["f_min", "f_max", "df_dt_max", "df_sign", "dt_init", "dt_ave", "tol", "with_kill"];
const DERIVED_KEYS = { air_anomaly: ["factor", "region"], obm_north: ["anomaly"], obm_south: ["anomaly"] };
const RECORD_KEYS = ["series_file", "series_var", "series_time_name"];

function derivedEq(name, p) {
  const a = sym(name) + "(\\tau)";
  const from = String(p.from || "");
  if (from === "air_anomaly") return `${a} = ${fmtNum(num(p.factor))}\\,\\overline{\\big(T^{\\mathrm{sl}}_{\\mathrm{ann}} - T^{\\mathrm{sl}}_{\\mathrm{ann},0}\\big)}^{\\,${p.region === "all" ? "x" : sub(p.region)}}`;
  if (from === "obm_north" || from === "obm_south") { const b = from === "obm_north" ? "n" : "s"; return `${a} = \\theta_{${b}}` + (p.anomaly ? ` - \\theta_{${b},\\mathrm{init}}` : ""); }
  return `${a} = \\text{?}\\qquad(\\texttt{from} = ${sub(from || "\"\"")})`;
}
// The record F of a coefficient in TeX (the file name, or F when none yet).
const recordSym = (p) => p.series_file ? `\\texttt{${String(p.series_file).replace(/^input\//, "").replace(/_/g, "\\_")}}${p.series_var ? ":" + sub(p.series_var) : ""}` : "F";
function coefEq(name, p) {
  const a = sym(name) + "(\\tau)", method = String(p.method), kind = kindOf(p);
  if (kind === "derived") return derivedEq(name, p) + `\\qquad\\text{(derived: set by the model each step)}`;
  if (kind === "constant") return `${a} = ${fmtNum(num(num(p.df_sign) >= 0 ? p.f_min : p.f_max))}`;
  if (kind === "record") return (p.norm ? `${a} = \\frac{F(\\tau) - \\bar F[t_0]}{\\bar F[t_1] - \\bar F[t_0]}` : `${a} = F(\\tau)`) + `,\\quad F = ${recordSym(p)}`;
  if (kind === "controller") return `\\dot ${sym(name)} = \\text{${method}}(\\dot V_{\\mathrm{ice}}),\\quad |\\dot ${sym(name)}| \\le ${fmtNum(p.df_dt_max)},\\quad ${sym(name)} \\in [${fmtNum(p.f_min)},\\,${fmtNum(p.f_max)}]`;
  return METHOD_EQ[method] ? METHOD_EQ[method](p, a) : `${a} = \\text{?}`;
}
// change the kind of a coefficient: the first method of the kind, a constant with f_min = f_max
function setKind(g, kind) {
  const p = coefPar(g.slice(COEF_PREFIX.length));
  const kv = { method: KINDS[kind][0] };
  if (kind === "constant") { kv.f_max = p.f_min; kv.df_dt_max = 0; }
  if (kind === "ramp" && num(p.df_dt_max) <= 0 && num(p.f_max) <= num(p.f_min)) { kv.f_max = num(p.f_min) + 1; kv.df_dt_max = 1e-3; }
  setMany(g, kv);
}

// The card of a coefficient in the equation: its kind and the parameters of the
// equation; a record also where it comes from (the file, the normalisation).
function coefCard(name, tl) {
  const g = COEF_PREFIX + name;
  if (!grp(g)) return h("div", { class: "fb-card" }, warn(`[${g}]: no such group and no coefficient template`), removeBtn(COEF_LIST_KEY, name));
  const { p, ser, f, kind } = coefOnRun(name, tl);
  const method = String(p.method);
  const kindSel = h("select", { onchange: (e) => setKind(g, e.target.value) }, ...Object.keys(KINDS).map(k => h("option", { value: k }, KIND_LABEL[k])));
  kindSel.value = kind;
  const head = h("div", { class: "fb-cardhead" }, texi(sym(name) + "(\\tau)"), h("code", {}, `[${g}]`), srcBadge(g), kindSel, removeBtn(COEF_LIST_KEY, name));
  const content = [head, newHint(g)];
  let shown = ["method", "units"];
  if (kind === "derived") {
    const from = String(p.from || "");
    const dk = DERIVED_KEYS[from] || ["factor", "region", "anomaly"];
    content.push(fields(g, ["from", ...dk, "units"]));
    shown.push("from", ...dk);
    content.push(from === "air_anomaly" ? note(`Set by the climate each update, after the atmosphere is composed: ${fmtNum(num(p.factor))} × the mean sea-level air temperature anomaly over ${p.region === "all" ? "the climate grid" : `region ${p.region}`}. Terms on it touch the ocean only.`)
      : from.startsWith("obm_") ? bipolarNote(`${from}: the ${from === "obm_north" ? "northern" : "southern"} box temperature of the ocean box model${p.anomaly ? " minus its initial value" : ""}, set each step`)
      : note("Choose what the coefficient is derived from."));
  } else if (kind === "constant") {
    const r = grp(g), v = num(num(p.df_sign) >= 0 ? p.f_min : p.f_max);
    const input = h("input", { type: "text", class: "fb-num", value: fmtNum(v), onchange: (e) => { const x = num(e.target.value); if (Number.isFinite(x)) setMany(g, { f_min: x, f_max: x }); else render(); } });
    content.push(h("div", { class: "fb-fields" }, fieldRow(g, "f_min", input, { label: "value", help: "The value of the coefficient at every forcing time (written as f_min = f_max)." }), field(g, "units")));
    shown.push("f_min", "f_max", "df_dt_max");
    if (r.given.has("f_max") && num(p.f_max) !== num(p.f_min) && num(p.df_dt_max) <= 0) content.push(note(`f_min = ${fmtNum(p.f_min)} and f_max = ${fmtNum(p.f_max)} differ; with no rate the coefficient stays at ${fmtNum(v)}.`));
  } else if (kind === "record") {
    content.push(fields(g, ["units"]));
    shown.push("norm_t0", "norm_t1", ...RECORD_KEYS);
  } else if (kind === "controller") {
    content.push(fields(g, ["method", ...CONTROLLER_KEYS, "units"], { method: { options: KINDS.controller } }));
    shown.push(...CONTROLLER_KEYS);
    content.push(note(`${method}: driven by the ice-volume response, no closed form. The curve is the fastest path, from ${num(p.df_sign) >= 0 ? "f_min" : "f_max"} at df_dt_max after dt_init.`));
  } else {
    const mk = METHOD_KEYS[method] || METHOD_KEYS["ramp-slope"];
    content.push(fields(g, ["method", ...mk, "units"], { method: { options: method === "const" ? ["const", ...KINDS.ramp] : KINDS.ramp } }));
    shown.push(...mk);
    const s = num(p.df_sign) >= 0 ? 1 : -1;
    const v = [`τ_e = max(0, τ − τ(t_init) − dt_init); s = ${s}`];
    if ((method === "const" || method === "ramp-slope") && num(p.df_dt_max) > 0) v.push(`from ${fmtNum(s > 0 ? p.f_min : p.f_max)} to ${fmtNum(s > 0 ? p.f_max : p.f_min)} over ${fmtNum(Math.abs(num(p.f_max) - num(p.f_min)) / num(p.df_dt_max))} yr`);
    if (tl.method === "const") v.push("the forcing time is constant: τ_e = 0, the coefficient stays at its start (tstep_method = rel advances it)");
    content.push(h("ul", { class: "fb-list" }, ...v.map(x => h("li", {}, x))));
  }
  if (num(p.win_t1) > num(p.win_t0)) content.push(vals(`window ${fmtNum(p.win_t0)} < τ < ${fmtNum(p.win_t1)}: ${name} → ${fmtNum(p.win_ref)} + ${fmtNum(p.win_fac)}·(${name} − ${fmtNum(p.win_ref)})`));
  if (num(p.sigma) > 0) content.push(vals(`noise: + ε, ε ~ N(0, ${fmtNum(p.sigma)}²) every step`));
  content.push(alsoSets(g, [...shown, "win_t0", "win_t1", "win_taper", "win_fac", "win_ref", "win_above", "sigma"]));
  if (kind === "record") content.push(recordFrom(name, p, ser, tl));
  else if (f && f.some(Number.isFinite)) {
    const rules = (kind === "ramp" || kind === "controller") ? [{ y: num(p.f_min), label: "f_min" }, { y: num(p.f_max), label: "f_max" }] : [];
    content.push(smallPlot([{ label: name, t: tl.t, f }], { y: name, rules }));
  }
  return h("div", { class: "fb-card" }, ...content);
}
// Where a piece of the equation gets its data, below its definition.
const fromBlock = (...content) => h("div", { class: "fb-from" }, h("div", { class: "fb-from-label" }, "from"), ...content);
// The record F of a coefficient: the file (and the variable of a netCDF file),
// the normalisation on it, the record over its full extent with the run's
// τ range as a yellow band and t₀/t₁ marked.
function recordFrom(name, p, ser, tl) {
  const g = COEF_PREFIX + name, r = grp(g);
  const netcdf = !!(p.series_file && D.series[p.series_file] && D.series[p.series_file].kind === "netcdf");
  const opts = { series_file: { options: ["", ...Object.keys(D.series).sort()] } };
  if (netcdf) opts.series_var = { options: Object.keys(D.series[p.series_file].vars) };
  const keys = ["series_file", ...(netcdf ? ["series_var", "series_time_name"] : [])];
  const normGiven = r.given.has("norm_t0") && r.given.has("norm_t1");
  const fixed = sourceDef(g) && sourceDef(g).keys.norm_t0 && sourceDef(g).keys.norm_t1;
  const cb = h("input", { type: "checkbox", disabled: !!fixed, title: fixed ? "the par file sets norm_t0 and norm_t1; -p cannot remove a key" : "",
    onchange: (e) => { if (e.target.checked) setMany(g, { norm_t0: val(g, "norm_t0", 0), norm_t1: num(val(g, "norm_t1", 0)) === num(val(g, "norm_t0", 0)) ? -21000 : val(g, "norm_t1") }); else { unset(g, "norm_t0"); unset(g, "norm_t1"); } } });
  cb.checked = normGiven;
  const content = [h("div", { class: "fb-fields" }, ...keys.map(k => field(g, k, opts[k] || {})),
    fieldRow(g, "norm", cb, { label: "normalise", help: "Rescale the record to 0 at t₀ and 1 at t₁ (a time, or two: the mean over the window), so the term adds the full anomaly at t₁." }),
    ...(normGiven ? [field(g, "norm_t0", { label: "t₀ (0)" }), field(g, "norm_t1", { label: "t₁ (1)" })] : []))];
  if (!p.series_file) content.push(warn("no record: the coefficient has no file"));
  else if (!ser) content.push(warn(`${p.series_file}${p.series_var ? ":" + p.series_var : ""}: not among the files exported to forcing.json (netCDF: set series_var)`));
  else {
    const vr = tl.method === "const" ? [{ x: tl.tc, label: "τ" }] : [];
    if (p.norm) { for (const t of asList(p.norm_t0).map(num)) vr.push({ x: t, label: "t0", stroke: "#393" }); for (const t of asList(p.norm_t1).map(num)) vr.push({ x: t, label: "t1", stroke: "#393" }); }
    if (p.norm) { const [v0, v1] = ser.norm; content.push(vals(`normalised: F̄[t₀] = ${fmtV(v0)}, F̄[t₁] = ${fmtV(v1)}`)); if (!Number.isFinite(v1 - v0) || Math.abs(v1 - v0) < 1e-12) content.push(warn("F̄[t₁] = F̄[t₀]: the normalisation divides by zero")); }
    const lo = ser.time[0], hi = ser.time[ser.time.length - 1];
    if (tl.tfmin < Math.min(lo, hi) || tl.tfmax > Math.max(lo, hi)) content.push(warn(`the record covers τ ∈ [${fmtNum(Math.min(lo, hi))}, ${fmtNum(Math.max(lo, hi))}] yr; the run reads τ ∈ [${fmtNum(tl.tfmin)}, ${fmtNum(tl.tfmax)}] yr (held at the end values outside)`));
    content.push(smallPlot([{ label: "F", t: ser.time, f: ser.f }], { x: "forcing time τ [yr]", y: p.norm ? "normalised" : "F", vrules: vr, band: tl.method !== "const" ? { x0: tl.tfmin, x1: tl.tfmax } : null }));
    content.push(vals(`${ser.n} points` + (ser.header ? `; ${ser.header.slice(0, 120)}` : "") + (tl.method !== "const" ? "; the run in yellow" : "")));
  }
  return fromBlock(...content);
}

// One term [term_<name>] as a record.
function termOf(name) {
  const g = "term_" + name;
  return { name, g, coef: String(val(g, TERM_COEF_KEY, "")), pattern: String(val(g, "pattern", "snapshots")), snap0: String(val(g, "snap0", "")), snap1: String(val(g, "snap1", "")),
    fields: asList(val(g, "fields", "")), gain: { tas: num(val(g, "f_ta", 1)), to: num(val(g, "f_to", 1)), so: num(val(g, "f_so", 1)) } };
}
const runTerms = () => termNames().filter(n => grp("term_" + n)).map(termOf);
// The snapshots of the run: the reference and the snap0/snap1 of the terms.
function runSnaps() {
  const snaps = [reference()];
  for (const t of runTerms()) if (t.pattern === "snapshots") for (const s of [t.snap0, t.snap1]) if (s && !snaps.includes(s)) snaps.push(s);
  return snaps;
}
// a snapshot group of the run: named, or created new by the equation
const snapGroup = (s) => s ? grp("snapshot_" + s) : null;
// the fields a snapshot sets
const setFields = (snap) => snapGroup(snap) ? SNAPFIELDS.filter(f => asList(val("snapshot_" + snap, f, "")).length) : [];
const climTdep = (b) => { const c = D.climdb[String(b)]; return !!(c && Array.isArray(c.time) && c.time[0] >= 1); };
const snapTdep = (snap) => SNAPFIELDS.some(f => asList(val("snapshot_" + snap, f, "")).some(climTdep));
const snapSym = (f, snap) => `${FSYM[f]}_{${sub(snap)}}` + (snapTdep(snap) ? "(\\tau)" : "");
// the coefficient of a term in TeX (warn colour when not defined yet); "" = 1
const coefTex = (t) => t.coef === "" ? "" : `${csym(t.coef)}(\\tau)`;
const wsym = (t) => t.coef === "" ? "" : coefTex(t) + "\\,";
const gsym = (x) => x === 1 ? "" : (x === -1 ? "-" : `${fmtNum(x)}\\,`);
function termRhs(t, f) {
  if (t.pattern === "replace") return `${FSYM[f]} = ${t.coef === "" ? "1" : coefTex(t)}\\quad\\text{at every depth, in place of the reference}`;
  return `${FSYM[f]} \\leftarrow ${FSYM[f]}` + termPart(t, f);
}
function termPart(t, f) {
  if (t.pattern === "uniform") {
    if (t.coef === "") return f === "pr" ? `\\,e^{\\beta_p\\,${fmtNum(t.gain.tas)}}` : ` ${t.gain[f] < 0 ? "-" : "+"} ${fmtNum(Math.abs(t.gain[f]))}`;
    if (f === "pr") return `\\,e^{\\beta_p\\,${t.gain.tas === 1 ? "" : "(" + fmtNum(t.gain.tas) + ")\\,"}${coefTex(t)}}`;
    const g = t.gain[f];
    return ` ${g < 0 ? "-" : "+"} ${gsym(Math.abs(g))}${coefTex(t)}`;
  }
  const s1 = t.snap1 || "?";
  if (f === "pr") return t.snap0 ? `\\,\\Big[1 + ${wsym(t)}\\Big(\\frac{${snapSym(f, s1)}}{${snapSym(f, t.snap0)}} - 1\\Big)\\Big]`
    : (t.coef === "" ? `\\,${snapSym(f, s1)}` : `\\,\\big[1 + ${wsym(t)}(${snapSym(f, s1)} - 1)\\big]`);
  return t.snap0 ? ` + ${wsym(t)}\\big[${snapSym(f, s1)} - ${snapSym(f, t.snap0)}\\big]` : ` + ${wsym(t)}${snapSym(f, s1)}`;
}
const PATTERN_HELP = { snapshots: "coefficient × (snap1 − snap0)", uniform: "coefficient × a gain, the same everywhere" };

// Structural warnings of a term (the equation); the snapshot blocks check the fields.
function termWarnings(t, tl) {
  const cl = climate();
  const sinfo = t.coef === "" ? null : (coefNames().includes(t.coef) ? coefOnRun(t.coef, tl) : null);
  const derived = !!(sinfo && sinfo.kind === "derived");
  const w = [];
  if (!t.fields.length) w.push("no fields: the term does nothing");
  if (derived && t.fields.some(f => f === "tas" || f === "pr")) w.push(`the derived coefficient "${t.coef}" is set after the atmosphere: this term may change to and so only`);
  if (t.pattern === "replace") {
    if (t.coef === "") w.push("replace needs a coefficient");
    if (t.fields.some(f => f === "tas" || f === "pr")) w.push("replace: to and so only");
  }
  if (t.pattern === "snapshots") {
    if (!t.snap1) w.push("no snap1: name the snapshot the term adds");
    if (derived) w.push("a derived coefficient cannot scale a snapshot pattern (the snapshots are composed before it is set)");
  }
  if (cl === "esm" && String(val("comps", "surface_method")) === "climate" && t.fields.some(f => f === "tas" || f === "pr")) w.push("the ESM supplies the smb (surface_method = climate): terms on tas/pr are refused");
  if (cl === "rembo" && t.fields.some(f => f === "tas" || f === "pr") && t.pattern !== "uniform") w.push("REMBO takes uniform tas terms only (its summer anomaly); snapshot terms on the atmosphere are ignored");
  return w;
}

// Define a coefficient a term names: the group of the par file or the
// library when there is one, else a new constant 1; appended to the list.
function defineCoef(n) {
  const g = COEF_PREFIX + n;
  const isnew = sourceOf(g) === "new" && !(S.g.has(g) && S.g.get(g).given.size);
  grp("forcing");
  const r = grp(g);
  if (isnew) { Object.assign(r.k, { method: "const", f_min: 1, f_max: 1 }); for (const k of ["method", "f_min", "f_max"]) r.given.add(k); }
  set("forcing", COEF_LIST_KEY, [...coefNames(), n]);
}
// The line of a term whose coefficient is not defined: what is missing, and
// the button that defines it.
function undefinedCoef(n) {
  const g = COEF_PREFIX + n, s = sourceOf(g);
  const there = s === "base" ? ` ([${g}] is in the par file)` : s === "library" ? ` ([${g}] is in the library)` : "";
  return h("div", { class: "fb-warn fb-define" }, `${usym(n)} is not defined yet${there}`,
    h("button", { class: "fb-btn", title: there ? `add [${g}] to the coefficients` : `a new [${g}]: the constant 1`, onclick: () => defineCoef(n) }, `define ${usym(n)}`));
}

function termCard(name, tl) {
  const g = "term_" + name;
  if (!grp(g)) return h("div", { class: "fb-card" }, warn(`[${g}]: no such group and no term template`), removeBtn("terms", name));
  const t = termOf(name);
  const content = [h("div", { class: "fb-cardhead" }, h("code", {}, `[${g}]`), srcBadge(g), h("span", { class: "fb-kind" }, t.pattern), h("span", { class: "fb-help" }, PATTERN_HELP[t.pattern] || ""), removeBtn("terms", name))];
  content.push(newHint(g), h("div", { class: "fb-fields" }, coefPicker(g, TERM_COEF_KEY, { label: "coefficient" }), field(g, "pattern")));
  if (t.pattern === "replace") content.push(bipolarNote("replace: the ocean field = the coefficient at every depth"));
  if (t.coef !== "" && !coefNames().includes(t.coef)) content.push(undefinedCoef(t.coef));
  const marks = {};
  if (t.pattern === "snapshots") {
    content.push(h("div", { class: "fb-fields" },
      snapPicker(g, "snap0", { label: "snap0 (coefficient = 0)", none: "none: snap1 is itself the anomaly or ratio" }),
      snapPicker(g, "snap1", { label: "snap1 (coefficient = 1)" })));
    for (const f of FIELDS) {
      const without = [t.snap0, t.snap1].filter(s => s && named("snapshot_" + s) && !setFields(s).includes(f));
      if (without.length) marks[f] = `not in snapshot ${without.join(", ")}`;
    }
  }
  if (t.pattern === "replace") for (const f of ["tas", "pr"]) marks[f] = "replace: to and so only";
  content.push(checkList(g, "fields", FIELDS, marks));
  if (t.pattern === "uniform") {
    const gk = ["tas", "to", "so"].filter(f => t.fields.includes(f)).map(f => GAIN[f]);
    if (gk.length) content.push(fields(g, gk));
    if (t.fields.includes("pr")) content.push(note(`pr follows the tas gain: × exp(β_p f_ta ${t.coef ? usym(t.coef) : "1"}), β_p = snapshots.f_p (fixed rule); leave pr out for no precipitation response.`));
  }
  // the term on its own: what it adds to (or puts in) each field
  const tf = FIELDS.filter(f => t.fields.includes(f));
  if (tf.length && (t.pattern !== "snapshots" || t.snap1))
    content.push(tex(tf.map(f => t.pattern === "replace" ? `${FSYM[f]} = ${t.coef === "" ? "1" : coefTex(t)}` : `${FSYM[f]}\\!:\\ ${f === "pr" ? "\\times" : ""}${termPart(t, f)}`).join(",\\qquad "), false));
  content.push(...termWarnings(t, tl).map(warn));
  if (t.pattern === "snapshots") {
    const snaps = [t.snap0, t.snap1].filter((s, i, a) => s && a.indexOf(s) === i);
    const role = (s) => `${s === t.snap0 && s === t.snap1 ? "snap0 and snap1" : s === t.snap0 ? "snap0" : "snap1"} of ${name}`;
    if (snaps.length) content.push(fromBlock(...snaps.map(s => snapshotBlock(s, role(s), tl))));
  }
  return h("div", { class: "fb-card" }, ...content);
}

// The fields a snapshot must set for the equation: the reference the
// atmosphere (tas, pr, zs) and every field of the terms on it; snap0/snap1
// the fields of their terms, plus zs when the pair is reduced to sea level.
function neededFields(snap) {
  const need = new Set();
  const terms = runTerms();
  if (snap === reference()) { if (climate() === "snapshots") { need.add("tas"); need.add("pr"); need.add("zs"); } for (const t of terms) for (const f of t.fields) if (FIELDS.includes(f)) need.add(f); }
  for (const t of terms) if (t.pattern === "snapshots" && (t.snap0 === snap || t.snap1 === snap)) {
    for (const f of t.fields) if (FIELDS.includes(f)) need.add(f);
    if (t.snap0 && (t.fields.includes("tas") || t.fields.includes("pr"))) need.add("zs");
  }
  return SNAPFIELDS.filter(f => need.has(f));
}
// a time-dependent field against the run's forcing time
function coverage(c, tl) {
  if (!c || !Array.isArray(c.time) || c.time[0] < 1) return null;
  const [, a, b] = c.time; const lo = Math.min(a, b), hi = Math.max(a, b);
  return { lo, hi, ok: tl.tfmin >= lo && tl.tfmax <= hi };
}
// The roles of a snapshot in the equation: the reference, snap0/snap1 of a term.
const snapRoles = (name) => [name === reference() ? "the reference" : null, ...runTerms().filter(t => t.pattern === "snapshots" && (t.snap0 === name || t.snap1 === name)).map(t => `${t.snap0 === name ? "snap0" : "snap1"} of ${t.name}`)].filter(Boolean);
// The entry of a group of input/clim_db.nml as a namelist block.
function climEntry(b, c) {
  const keys = ["filename", "name", "units", "time", ...(c.grid_src ? ["grid_src"] : [])];
  const w = Math.max(...keys.map(k => k.length));
  return [h("pre", { class: "fb-cfg" }, `&${b}\n${keys.map(k => `    ${k.padEnd(w)} = ${fmtNml(c[k])}`).join("\n")}\n/`),
    h("p", { class: "fb-note" }, "defined in input/clim_db.nml; the format: ", h("a", { href: "database-namelists.html" }, "database namelists"))];
}
// The clim_db groups of a snapshot as a table: field, group, variable, file,
// τ range; a click on the group shows its entry in a row below (kept open
// across redraws, by id).
function climTable(id, rows) {
  const body = [];
  for (const { f, b, c, cov } of rows) {
    const key = `clim-${id}-${f}-${b}`;
    const det = h("tr", { class: "fb-tbl-det", hidden: !S.open.has(key) },
      h("td", { colspan: 5 }, c ? climEntry(b, c) : note(`"${b}" is not a group of input/clim_db.nml`)));
    const toggle = () => { det.hidden = !det.hidden; if (det.hidden) S.open.delete(key); else S.open.add(key); };
    body.push(h("tr", {}, h("td", {}, f), h("td", {}, h("button", { class: "fb-link", title: "show the entry in input/clim_db.nml", onclick: toggle }, b)),
      h("td", {}, c ? String(c.name) : "–"), h("td", { class: "fb-tbl-file" }, c ? c.filename.replace(/^fesmdata\//, "") : "not in input/clim_db.nml"),
      h("td", {}, !c ? "" : cov ? `${fmtNum(cov.lo)} to ${fmtNum(cov.hi)} yr` : "static")), det);
  }
  return h("table", { class: "fb-tbl" }, h("thead", {}, h("tr", {}, ...["field", "group", "variable", "file", "τ range"].map(x => h("th", {}, x)))), h("tbody", {}, ...body));
}
// The block of [snapshot_<name>] inside the card that names it (role: what it
// is there): the clim_db groups of each field as chips, marked where the
// equation needs it; a table of where each group points (variable, file,
// τ range) with its clim_db entry on a click; how the fields are read.
// The same snapshot named twice is the same group: edits apply to both.
function snapshotBlock(name, role, tl) {
  const g = "snapshot_" + name;
  const roles = snapRoles(name);
  const head = (...extra) => h("div", { class: "fb-cardhead" }, h("code", {}, `[${g}]`), srcBadge(g), h("span", { class: "fb-help" }, role, ...extra));
  if (!grp(g)) return h("div", { class: "fb-sub-card" }, head(), warn(`[${g}]: no group and no snapshot template`));
  const need = neededFields(name);
  const rows = [], tbl = [], warns = [];
  for (const f of SNAPFIELDS) {
    const groups = asList(val(g, f, ""));
    const needed = need.includes(f);
    if (!groups.length && !needed) continue;
    for (const b of groups) {
      const c = D.climdb[b] || null, cov = c ? coverage(c, tl) : null;
      if (!c) warns.push(`${f}: "${b}" is not a group of input/clim_db.nml`);
      if (cov && !cov.ok) warns.push(`${f}: ${b} covers τ ∈ [${fmtNum(cov.lo)}, ${fmtNum(cov.hi)}] yr; the run reads τ ∈ [${fmtNum(tl.tfmin)}, ${fmtNum(tl.tfmax)}] yr (held at the end slices outside)`);
      tbl.push({ f, b, c, cov });
    }
    const label = !needed ? h("span", {}, h("code", {}, f), h("span", { class: "fb-help" }, " (not needed)"))
      : groups.length ? h("code", {}, f + " ✓") : h("code", { class: "fb-miss" }, f + ": missing");
    rows.push(chipList(g, f, climdbNames(), { label, help: "" }));
    if (needed && !groups.length) warns.push(`${f} is needed by ${roles.join(", ")} but not set`);
  }
  const placeholders = tbl.some(x => x.c && /\{\w+\}/.test(x.c.filename));
  const uniq = [...new Set(warns)];
  return h("div", { class: "fb-sub-card" },
    head(snapTdep(name) ? "; depends on τ" : ""), newHint(g),
    h("div", { class: "fb-fields" }, ...rows),
    rows.length ? note("Each field names its clim_db group(s): tas one monthly group, or two for annual + summer; pr one, or two summed (snowfall + rainfall); the others one.") : null,
    rows.length ? climdbHint() : null,
    tbl.length ? climTable(`${name}-${role}`, tbl) : note("No fields set: add a clim_db group per field, or pick a snapshot of the library."),
    placeholders ? note("In the file names, {domain} and {grid_name} are filled from the run's domain and grid, {grid_global} from [snapshots] grid_global, {exp} from the snapshot's exp.") : null,
    details(`snap-${name}-${role}`, "how the fields are read (exp, monthly, smoothing)", fields(g, ["exp", "monthly", "clim_sigma", "ocn_sigma"]), alsoSets(g, [...SNAPFIELDS, "exp", "monthly", "clim_sigma", "ocn_sigma"])),
    ...uniq.map(warn));
}

// The explainer at the top of the page: what a snapshot, a term and a
// coefficient are (left); one example with two terms, its pieces named, on
// both equations of the atmosphere (right). KaTeX needs literal colours: the
// page's muted, accent, base green and warn.
const GUIDE_C = { ref: "#6b6861", coef: "#1f5f8b", pat: "#2a5d24", par: "#6250d6", term: "#9a4b00" };
function guideEquation() {
  const c = (k, x) => `\\textcolor{${GUIDE_C[k]}}{${x}}`;
  const ub = (k, x, lab) => `\\underbrace{${c(k, x)}}_{${c(k, lab)}}`;
  const ob = (x, lab = "\\text{term}") => `\\overbrace{${x}}^{${c("term", lab)}}`;
  const T = (n) => `T^{\\mathrm{sl}}_{\\mathrm{${n}}}`, P = (n) => `P^{\\mathrm{sl}}_{\\mathrm{${n}}}`;
  const alpha = ub("coef", "\\alpha(\\tau)", "\\text{coefficient}");
  // the uniform term: the gain f (and β_p) are parameters, β(τ) the coefficient,
  // the pattern the unit field (named on the term's brace)
  const beta = ub("coef", "\\beta(\\tau)", "\\text{coefficient}");
  const uni = "\\text{term (pattern: uniform)}";
  // one line each: the explainer gives the equations two thirds of the width
  const tas = `T^{\\mathrm{sl}}(x,m,\\tau) &= ${ub("ref", T("pd"), "\\text{reference}")}`
    + ` + ${ob(`${alpha}\\,${ub("pat", `\\big[${T("lgm")} - ${T("picontrol")}\\big]`, "\\text{pattern: snapshots}")}`)}`
    + ` + ${ob(`${ub("par", "f", "\\text{parameter}")}\\;${beta}`, uni)}`;
  const pr = `P^{\\mathrm{sl}}(x,m,\\tau) &= ${ub("ref", P("pd"), "\\text{reference}")}`
    + ` \\cdot ${ob(`\\Big[1 + ${alpha}\\,${ub("pat", `\\Big(\\frac{${P("lgm")}}{${P("picontrol")}} - 1\\Big)`, "\\text{pattern: snapshots (ratio)}")}\\Big]`)}`
    + ` \\cdot ${ob(`\\exp\\!\\big(${ub("par", "\\beta_p\\, f", "\\text{parameters}")}\\;${beta}\\big)`, uni)}`;
  return `\\begin{aligned} ${tas} \\\\[1.4em] ${pr} \\end{aligned}`;
}
function explainer() {
  return h("div", { class: "fb-explainer" },
    h("p", { class: "fb-note fb-primer" },
      "A snapshot is a named set of climate fields (tas, pr, zs, to, so), read from files. "
      + "A term is always coefficient × pattern: the coefficient c(τ) is one number at each forcing time (a constant, a record, a ramp, ...), the pattern a field that carries the spatial structure and the units. "
      + "A parameter is a fixed number set in the par file (the gain f of a uniform term, β_p, the lapse rates); it does not change with τ. "
      + "In a snapshots pattern the field is S₁ − S₀ in the units of the field, and c has no units (0 at S₀, 1 at S₁); in a uniform pattern the field is the unit field times the gain f, f in units of the field per unit of c, so an offset of c kelvin is c in K times the gain 1. "
      + "The terms add to the reference (pr multiplies)."
      + (climate() === "esm" ? " With the ESM backend, the climate of the ESM takes the place of the reference." : "")),
    h("div", { class: "fb-guide" }, tex(guideEquation()),
      note("An example with two terms: a glacial-index term as in the Greenland paleo setup and a uniform term as in the Pyrenees or hysteresis setups. One term sets both equations: tas adds, pr multiplies, by the fixed rule.")));
}
function secEquation(tl) {
  const cl = climate();
  const content = [];
  content.push(badgeLegend());
  if (cl === "esm") {
    content.push(card("[esm]", null, h("span", { class: "fb-help" }, "the ESM supplies the climate; the terms add to it"),
      fields("esm", ["experiment", "esm_name", "use_esm", "use_hist", "use_proj"]),
      note("The ESM backend composes the atmosphere and the ocean from its own anomalies (a reference climatology plus the ESM anomaly over its reference period). No reference snapshot and no snapshot terms here: uniform terms on tas/pr unless the ESM supplies the smb, ocean terms at the shelf base. The periods and the ESM files: [esm] and the ESM par file (climate-esm.html)."),
      h("p", { class: "fb-note" }, h("a", { href: "climate-esm.html" }, "the ESM backend"))));
  } else {
    content.push(h("div", { class: "fb-card fb-card-flat" },
      h("div", { class: "fb-fields" }, snapPicker("forcing", "reference")),
      note(cl === "snapshots" ? "tas, to, so add; pr multiplies (fixed by field). Applied in this order: the snapshot terms, the uniform terms on tas/pr, the derived coefficients, the ocean terms."
        : "REMBO supplies the atmosphere and takes the uniform tas terms as its summer anomaly (ΔT_sum = Σ f_ta α, ΔT_ann = 1.3 ΔT_sum); the snapshot backend supplies the ocean from the reference and the ocean terms."),
      fromBlock(snapshotBlock(reference(), "the reference", tl),
        fields("snapshots", ["lapse", "f_p"], { lapse: { label: "lapse (Γ)" }, f_p: { label: "f_p (β_p)" } }),
        note("Γ (annual, summer) and β_p: the sea-level reduction of every snapshot, the uniform terms and the elevation transform use the same values."))));
  }
  const cnames = coefNames();
  content.push(h("h3", { class: "fb-sub" }, "Terms ", h("span", { class: "fb-hint" }, "coefficient × pattern on explicit fields, added to the reference")));
  content.push(...termNames().map(n => termCard(n, tl)));
  content.push(addBox("terms", "term_", "term", (g) => { const r = grp(g); Object.assign(r.k, { snap0: "picontrol", snap1: "lgm", fields: ["tas", "pr"], [TERM_COEF_KEY]: cnames[0] || "" }); for (const k of [TERM_COEF_KEY, "snap0", "snap1", "fields"]) r.given.add(k); }));
  content.push(h("h3", { class: "fb-sub" }, "Coefficients ", h("span", { class: "fb-hint" }, "α(τ), one number at each forcing time: a constant, a record, a ramp, a controller, or derived from the model")));
  content.push(...cnames.map(n => coefCard(n, tl)));
  content.push(addBox(COEF_LIST_KEY, COEF_PREFIX, "coefficient", (g) => { const r = grp(g); r.k.method = "series"; r.given.add("method"); }));
  return section("equation", "The equation", "reference + Σ coefficient × pattern, each piece with its source", ...content);
}

// ---------------------------------------------- left: after the composition

function secAfter() {
  const content = [];
  content.push(
    card("[taper]", null, h("span", { class: "fb-help" }, "the atmosphere towards the reference with a weight field"), fields("taper", ["active", ...(val("taper", "active", false) ? ["weight", "outside"] : [])])),
    card("[bias]", null, h("span", { class: "fb-help" }, "an additive temperature field and a precipitation factor field"), fields("bias", ["active", ...(val("bias", "active", false) ? ["dtas", "fpr", "f_tas", "f_pr"] : [])], { dtas: { chips: climdbNames(), single: true }, fpr: { chips: climdbNames(), single: true } }),
      val("bias", "active", false) ? climdbHint() : null));
  const om = String(val("orbit", "method", "fixed_time"));
  content.push(card("[orbit]", null, h("span", { class: "fb-help" }, "insolation of the surface models"), fields("orbit", ["method", ...(om === "fixed_time" ? ["time"] : om === "fixed_orbit" ? ["ecc", "obl", "per"] : [])])));
  const m = String(val("barysealevel", "method", "fastiso"));
  const sl = [fields("barysealevel", ["method", ...(m === "const" ? ["bsl_init"] : m === "file" ? ["sl_path", "sl_name"] : [])], { sl_path: { options: ["none", ...Object.keys(D.series).sort()] } })];
  if (m === "file") {
    const path = String(val("barysealevel", "sl_path", "")), vname = String(val("barysealevel", "sl_name", ""));
    if (!seriesOf(path, vname)) sl.push(warn(`${path}${vname && vname !== "none" ? ":" + vname : ""}: not among the files exported to forcing.json (netCDF: set sl_name)`));
  }
  content.push(card("[barysealevel]", null, h("span", { class: "fb-help" }, "sea level"), ...sl));
  if (inRun("ice_bnd") && sourceOf("ice_bnd")) {
    const im = String(val("ice_bnd", "method")), on = val("ice_bnd", "active", false);
    const ib = [fields("ice_bnd", ["active", ...(on ? ["method", "H_ice", ...(im === "scaled" ? ["coef", "slice_time"] : [])] : [])], { coef: { options: ["", ...coefNames()] }, H_ice: { chips: climdbNames(), single: true } })];
    if (on) ib.push(climdbHint());
    if (on && im === "scaled" && !coefNames().includes(String(val("ice_bnd", "coef", "")))) ib.push(warn(`ice_bnd.coef = "${val("ice_bnd", "coef", "")}" is not a coefficient of the equation`));
    content.push(card("[ice_bnd]", null, h("span", { class: "fb-help" }, "prescribed ice in the fixed and relaxation cells"), ...ib));
  }
  return section("after", "After the composition", "to the surface (always); transforms; also at the forcing time", ...content);
}

// ------------------------------------------------- right: the equation chain

// The evaluation time: one slider and its readout, created once and kept
// across renders, so that a drag survives the redraw of the blocks below it.
let SLIDER = null, TAUVAL = null;
const tauText = (tl, t) => `t = ${fmtNum(+t.toPrecision(6))} yr → τ = ${fmtNum(+tauAt(tl, t).toPrecision(6))} yr`;
function timeSlider(tl) {
  if (!SLIDER) {
    SLIDER = h("input", { type: "range" });
    TAUVAL = h("span", { class: "fb-tauval" });
    SLIDER.addEventListener("input", () => { S.tsel = Number(SLIDER.value); renderBelow(); });
  }
  const lo = Math.min(tl.t0, tl.t1), hi = Math.max(tl.t0, tl.t1), t = evalTime(tl);
  Object.assign(SLIDER, { min: lo, max: hi, step: (hi - lo) / 400 || 1 });
  SLIDER.value = t;
  TAUVAL.textContent = tauText(tl, t);
  return h("div", { class: "fb-taurow" }, SLIDER, TAUVAL);
}
// Block 0: the forcing time and the slider of the evaluation time.
function rightTime(tl) {
  const timeEq = tl.method === "const" ? `\\tau = t_{\\mathrm{const}} = ${fmtK(tl.tc)}` : tl.method === "rel" ? `\\tau = t` : `\\tau = t_{\\mathrm{cal}} - ${fmtNum(tl.cref)}`;
  return eqblock("0 · Forcing time",
    tex(`${timeEq},\\qquad t \\in [${fmtK(tl.t0)},\\ ${fmtK(tl.t1)}]\\ \\mathrm{yr}\\quad(\\texttt{${tl.method}})`),
    timeSlider(tl),
    vals("Every coefficient, time-dependent snapshot, sea-level file and the orbit (forcing_time) read τ. The values below are at this τ."));
}
// Blocks 1 to 6: the equation chain at the evaluation time.
function rightEquations(tl) {
  const t = evalTime(tl), tau = tauAt(tl, t);
  const out = [];
  const cl = climate(), snames = coefNames(), ref = reference();
  const terms = runTerms();
  const sv = {};
  for (const n of snames) { const s = coefOnRun(n, tl); sv[n] = !s.p || s.kind === "derived" ? null : atT(tl, s.f, t); }
  const tauS = `${fmtNum(+tau.toPrecision(6))} yr`;

  // the coefficients the terms name but the equation does not define
  const undef = [...new Set(terms.map(x => x.coef).filter(n => n !== "" && !snames.includes(n)))];
  if (snames.length || undef.length) {
    const b = eqblock("1 · Coefficients");
    const used = new Set(terms.map(x => x.coef));
    const ibCoef = inRun("ice_bnd") && val("ice_bnd", "active", false) && String(val("ice_bnd", "method")) === "scaled" ? String(val("ice_bnd", "coef", "")) : "";
    for (const n of snames) {
      const p = coefPar(n);
      if (!p) { b.append(warn(`coefficient "${n}": no group [${COEF_PREFIX}${n}]`)); continue; }
      const at = sv[n] === null || !Number.isFinite(sv[n]) ? "" : `\\qquad ${sym(n)}(${fmtK(+tau.toPrecision(6))}) = ${fmtV(sv[n])}`;
      const unused = used.has(n) ? "" : `\\qquad\\text{(used by no term${n === ibCoef ? "; scales ice\\_bnd" : ""})}`;
      b.append(tex(coefEq(n, p) + at + unused));
    }
    for (const n of undef) b.append(tex(`${csym(n)}(\\tau)\\qquad{\\color{${WARN_TEX}}\\text{not defined yet}}`));
    out.push(b);
  }

  if (withSnaps()) {
    const lapse = asList(val("snapshots", "lapse", [0.008, 0.0065]));
    out.push(eqblock("2 · Snapshots, reduced to sea level once",
      tex(`T^{\\mathrm{sl}}_j = T_j + \\Gamma_m z_j,\\qquad P^{\\mathrm{sl}}_j = P_j\\,e^{\\beta_p \\Gamma_m z_j},\\qquad \\Theta_j(z),\\ S_j(z)\\ \\text{as given}`),
      vals(`Γ = ${fmtNum(num(lapse[0]))}, ${fmtNum(num(lapse[1]))} K/m (annual, summer; monthly cosine between them), β_p = ${fmtNum(num(val("snapshots", "f_p", 0.05)))} /K; j ∈ {${runSnaps().map(s => s + (snapTdep(s) ? "(τ)" : "")).join(", ")}}; reference S₀ = ${ref}`)));
    const b3 = eqblock(cl === "rembo" ? "3 · Composition (REMBO: the atmosphere; the snapshots: the ocean)" : "3 · Composition at sea level: reference + Σ coefficient × pattern");
    if (!snapGroup(ref)) b3.append(warn(`the reference snapshot "${ref}" has no group [snapshot_${ref}]`));
    else {
      const refF = setFields(ref);
      for (const f of FIELDS) {
        const tf = terms.filter(t => t.fields.includes(f));
        if (!refF.includes(f) && !tf.length) continue;
        if (cl === "rembo" && (f === "tas" || f === "pr")) continue;
        const rp = tf.find(t => t.pattern === "replace");
        let e = `${FSYM[f]}(${f === "tas" || f === "pr" ? "x,m" : "x,z"},\\tau) = ` + (rp ? (coefTex(rp) || "1") : snapSym(f, ref));
        for (const t of tf) if (t.pattern !== "replace") e += termPart(t, f);
        b3.append(tex(e));
      }
      if (cl === "rembo") { const u = terms.filter(t => t.pattern === "uniform" && t.fields.includes("tas")); if (u.length) b3.append(tex(`\\Delta T_{\\mathrm{sum}} = ${u.map(t => gsym(t.gain.tas) + (coefTex(t) || "1")).join(" + ")},\\qquad \\Delta T_{\\mathrm{ann}} = 1.3\\,\\Delta T_{\\mathrm{sum}}`)); }
      for (const f of ["to", "so"]) if (terms.filter(t => t.pattern === "replace" && t.fields.includes(f)).length > 1) b3.append(warn(`more than one replace term on ${f}`));
    }
    const num3 = terms.map(t => t.coef === "" ? `${t.name}: coefficient 1` : !snames.includes(t.coef) ? `${t.name}: ${t.coef} not defined yet` : sv[t.coef] === null ? `${t.name}: ${t.coef} set at run time` : `${t.name}: ${t.coef} = ${fmtV(sv[t.coef])}` + (t.pattern === "uniform" && t.fields.includes("tas") ? ` → ΔT = ${fmtV(t.gain.tas * sv[t.coef])} K` : ""));
    if (num3.length) b3.append(vals(`at τ = ${tauS}:  ` + num3.join("  ·  ")));
    out.push(b3);
    if (cl === "snapshots") out.push(eqblock("4 · To the model surface z_s(x,t)", tex(`T(x,m) = T^{\\mathrm{sl}} - \\Gamma_m z_s,\\qquad P(x,m) = P^{\\mathrm{sl}}\\,e^{-\\beta_p \\Gamma_m z_s}`)));
  } else if (cl === "esm") {
    const useEsm = val("esm", "use_esm", false), exp = String(val("esm", "experiment", "ctrl"));
    const b = eqblock("2 · The ESM climate");
    if (useEsm && exp !== "ctrl") b.append(tex(`X(t) = X_{\\mathrm{ref}} + \\big(X_{\\mathrm{ESM}}(t) - \\overline{X_{\\mathrm{ESM}}}^{\\,\\mathrm{time\\_esm\\_ref}}\\big),\\quad \\text{historical to } ${fmtNum(asList(val(tl.T, "time_hist", [0, 0]))[1])}\\text{, then } \\texttt{${exp}}`),
      vals(`${val("esm", "esm_name", "")}: air temperature and smb anomalies are differences, precipitation a ratio, the ocean annual means at the shelf base.`));
    else if (!useEsm && exp === "1pctCO2") b.append(tex(`\\Delta T_a(t) = f_{\\mathrm{polar}}\\,\\min\\!\\big(17\\,\\mathrm{K}\\,\\tfrac{t - 1995}{305},\\ \\Delta T_{\\mathrm{thr}}\\big),\\qquad \\Delta\\Theta = f_{\\mathrm{ocn}}\\,\\Delta T_a`),
      vals(`f_polar = ${fmtNum(val("esm", "f_polar"))}, ΔT_thr = ${fmtNum(val("esm", "dT_threshold"))}, f_ocn = ${fmtNum(val("esm", "f_ocn"))}: homogeneous warming (use_esm = False, experiment 1pctCO2).`));
    else b.append(vals(`experiment = "${exp}", use_esm = ${useEsm}: the reference climatology only (spin-up).`));
    out.push(b);
    if (terms.length) { const b3 = eqblock("3 · Terms on the ESM result"); for (const t of terms) for (const f of t.fields) if (FIELDS.includes(f)) b3.append(tex(`\\text{${t.name}:}\\quad ` + termRhs(t, f))); out.push(b3); }
  }

  const tr = [];
  if (val("taper", "active", false)) tr.push(tex(`X \\leftarrow X_0(z_s) + w(x)\\,\\big(X - X_0(z_s)\\big)` + (String(val("taper", "outside")) === "mean" ? ` + (1-w)\\,\\overline{(X - X_0)}^{\\,w=1}` : "") + `\\qquad X \\in \\{T, P\\}`),
    vals(`taper: w = ${val("taper", "weight", "")} (1: the composed climate, 0: the reference at the current surface)`));
  if (val("bias", "active", false)) {
    const tdep = climTdep(val("bias", "dtas", "")) || climTdep(val("bias", "fpr", ""));
    const a = tdep ? "(x,m,\\tau)" : "(x,m)";
    tr.push(tex(`T \\leftarrow T + ${gsym(num(val("bias", "f_tas")))}b_T${a},\\qquad P \\leftarrow P\\,\\big(1 + ${gsym(num(val("bias", "f_pr")))}(b_P${a} - 1)\\big)`), vals(`bias: b_T = ${val("bias", "dtas", "")}, b_P = ${val("bias", "fpr", "")}${tdep ? " (a time axis: sliced at τ)" : ""}`));
  }
  out.push(eqblock("5 · Transforms", ...(tr.length ? tr : [vals("none active (taper, bias)")])));

  const b6 = eqblock("6 · Also at the forcing time");
  const om = String(val("orbit", "method", "fixed_time"));
  b6.append(om === "forcing_time" ? tex(`I = I\\big(\\varphi, d;\\ \\mathrm{orbit}(\\tau)\\big)\\quad\\text{Laskar (2004)}`)
    : om === "fixed_time" ? tex(`I = I\\big(\\varphi, d;\\ \\mathrm{orbit}(${fmtK(num(val("orbit", "time")))}\\ \\mathrm{yr})\\big)\\quad\\text{fixed, Laskar (2004)}`)
    : om === "fixed_orbit" ? tex(`I = I(\\varphi, d;\\ e, \\varepsilon, \\varpi),\\qquad e = ${fmtNum(val("orbit", "ecc"))},\\ \\varepsilon = ${fmtNum(val("orbit", "obl"))}^\\circ,\\ \\varpi = ${fmtNum(val("orbit", "per"))}^\\circ`)
    : warn(`orbit.method = "${om}": forcing_time, fixed_time or fixed_orbit`));
  const m = String(val("barysealevel", "method", "fastiso"));
  if (m === "const") b6.append(tex(`\\eta = ${fmtNum(val("barysealevel", "bsl_init"))}\\ \\mathrm{m}`));
  else if (m === "fastiso") b6.append(tex(`\\eta(t) = \\eta_0 - \\Delta V_{\\mathrm{ice}}(t)/A_{\\mathrm{ocean}}\\quad\\text{(the FastIsostasy domains)}`));
  else if (m === "file") { const ser = seriesOf(String(val("barysealevel", "sl_path", "")), String(val("barysealevel", "sl_name", ""))); b6.append(tex(`\\eta = \\eta_{\\mathrm{file}}(\\tau)` + (ser ? `\\qquad \\eta(${fmtK(+tau.toPrecision(6))}) = ${fmtV(interp(ser, tau))}\\ \\mathrm{m}` : ""))); }
  if (inRun("ice_bnd") && val("ice_bnd", "active", false)) {
    const im = String(val("ice_bnd", "method")), s = String(val("ice_bnd", "coef", ""));
    b6.append(im === "scaled" ? tex(`H_b(x) = ${s ? sym(s) + "(\\tau)" : "?"}\\,H(x, ${fmtK(num(val("ice_bnd", "slice_time")))})`) : tex(`H_b(x) = H(x, \\tau)`));
  }
  out.push(b6);
  return out;
}

// ------------------------------------------------------ right: over the run

function rightRun(tl) {
  const t = evalTime(tl);
  const terms = runTerms();
  const content = [];
  const ser = coefNames().map(n => ({ n, s: coefOnRun(n, tl) })).filter(x => x.s.f && x.s.f.some(Number.isFinite)).map(x => ({ label: x.n, t: tl.t, f: x.s.f }));
  if (ser.length) content.push(linePlot(ser, { y: "coefficients", tmark: t }));
  const uni = terms.filter(x => x.pattern === "uniform");
  if (uni.length) {
    const sums = {};
    for (const x of uni) {
      const w = termWeight(x, tl);
      if (!w) continue;
      for (const [f, lab] of [["tas", climate() === "rembo" ? "ΔT_sum [K]" : "ΔT_a [K]"], ["to", "ΔΘ [K]"], ["so", "ΔS [psu]"]]) if (x.fields.includes(f)) sums[lab] = (sums[lab] || tl.t.map(() => 0)).map((v, i) => v + x.gain[f] * w[i]);
    }
    const s = Object.entries(sums).map(([label, f]) => ({ label, t: tl.t, f }));
    if (s.length) content.push(linePlot(s, { y: "uniform anomalies", height: 160, tmark: t }));
  }
  if (String(val("barysealevel", "method")) === "file") {
    const sl = seriesOf(String(val("barysealevel", "sl_path", "")), String(val("barysealevel", "sl_name", "")));
    if (sl) content.push(linePlot([{ label: "η [m]", t: tl.t, f: tl.tf.map(x => interp(sl, x)) }], { y: "sea level [m]", height: 140, tmark: t }));
  }
  if (!content.length) content.push(vals(tl.method === "const" ? "The forcing time is constant: nothing changes over the run." : "No coefficient with a curve in the run (derived coefficients are set by the model)."));
  return eqblock("Over the run (model time; the line: the evaluation time)", ...content);
}

// ------------------------------------------------- right: the run written out

// The groups of the run, in par-file order: the time groups, the forcing,
// the named groups the equation references, what comes after.
function runGroups() {
  const T = tgroup(), order = [];
  const add = (g, always = true) => {
    if (!g || order.includes(g) || !sourceOf(g)) return;
    if (!always && sourceOf(g) === "new" && !(S.g.has(g) && [...S.g.get(g).given].some(k => changed(g, k)))) return;
    order.push(g);
  };
  add("ctrl"); if (T !== "ctrl") add(T);
  if (climate() === "esm") add("esm");
  add("forcing");
  if (withSnaps()) { add("snapshots"); for (const s of runSnaps()) add("snapshot_" + s); }
  for (const n of coefNames()) add(COEF_PREFIX + n);
  for (const n of termNames()) add("term_" + n);
  for (const g of ["taper", "bias", "orbit", "barysealevel", "ice_bnd"]) add(g, false);
  return order;
}
function fmtP(v) {
  if (typeof v === "boolean") return v ? "True" : "False";
  if (typeof v === "number") return fmtNum(v);
  const s = String(v);
  return /^[\w.\/{}:~&|,-]*$/.test(s) ? s : JSON.stringify(s);
}
function fmtNml(v) {
  if (Array.isArray(v)) return v.length ? v.map(fmtNml).join(" ") : "\"\"";
  if (typeof v === "boolean") return v ? "True" : "False";
  if (typeof v === "number") return fmtNum(v);
  return JSON.stringify(String(v));
}
// Tab "par file": the groups of the run as the model reads them; the keys
// each group sets, highlighted where the builder changed them.
function parFile() {
  const lines = [];
  const esc = (s) => s.replace(/&/g, "&amp;").replace(/</g, "&lt;");
  for (const g of runGroups()) {
    const r = grp(g), src = sourceOf(g);
    const keys = schemaKeys(g).filter(k => r.given.has(k));
    lines.push(`<span class="g">&amp;${g}</span>` + (src === "library" ? `   <span class="c">! from the forcing library (${esc(D.library.files.join(", "))})</span>` : src === "new" ? `   <span class="c">! new: add to the par file or to input/forcing_db.nml</span>` : ""));
    const w = Math.max(...keys.map(k => k.length), 4);
    for (const k of keys) {
      const ch = changed(g, k), ref = refOf(g, k);
      const c = ch ? (ref && ref.src !== "default" ? `! was ${fmtNml(ref.v)}` : "") : "";
      lines.push(`<span${ch ? ' class="ov"' : ""}>    ${k.padEnd(w)} = ${esc(fmtNml(r.k[k]))}${c ? `   <span class="c">${esc(c)}</span>` : ""}</span>`);
    }
    lines.push("/", "");
  }
  if (S.extra.length) lines.push(`<span class="c">! -p outside the forcing groups, passed through to the command:</span>`, ...S.extra.map(x => `<span class="c">!   ${esc(x)}</span>`));
  return lines.join("\n");
}
// Tab "-p": what differs from the base par file (and the library groups it names).
function overrides() {
  const out = [], blocks = [], problems = [];
  for (const g of runGroups()) {
    const r = grp(g), src = sourceOf(g);
    const ids = schemaKeys(g).filter(k => changed(g, k));
    if (!ids.length) continue;
    if (src === "new") {
      blocks.push(`&${g}\n${ids.map(k => `    ${k.padEnd(16)} = ${fmtNml(r.k[k])}`).join("\n")}\n/`);
      problems.push(`[${g}] is new: -p cannot create it. Add the block below to the par file or to input/forcing_db.nml.`);
      continue;
    }
    for (const k of ids) out.push(Array.isArray(r.k[k]) ? `${g}.${k}==${r.k[k].map(fmtP).join(",")}` : `${g}.${k}=${fmtP(r.k[k])}`);
    if (src === "library") problems.push(`[${g}] comes from the library: -p on it needs runme ≥ 0.5.26 (it copies the group into the staged par file).`);
  }
  return { out: [...out, ...S.extra], blocks, problems };
}
function rightOutput() {
  const tabs = [["par", "par file"], ["p", "-p and the runme command"]];
  const head = h("div", { class: "fb-tabs" }, ...tabs.map(([k, lab]) => h("button", { class: S.tab === k ? "on" : "", onclick: () => { S.tab = k; renderRight(); } }, lab)));
  const sel = h("select", { onchange: (e) => { S.base = e.target.value; render(); } }, ...Object.keys(D.bases).map(b => h("option", { value: b }, b)));
  sel.value = S.base;
  const body = [h("div", { class: "fb-fields" }, h("div", { class: "fb-field" }, h("label", { title: "runme -n" }, h("code", {}, "base par file (-n)")), h("span", { class: "fb-val" }, sel), h("span", { class: "fb-help" }, "The par file the run is written against: the -p list is what differs from it.")))];
  if (S.tab === "par") {
    body.push(note("The forcing groups of the run as the model reads them: the base par file, the library groups the equation names, new groups. Highlighted: what the builder changed."),
      h("pre", { class: "fb-cfg", html: parFile() }));
  } else {
    const { out, blocks, problems } = overrides();
    const outdir = S.outdir || `output/${S.preset || S.base.replace(/^yelmox_/, "").toLowerCase() + "-run"}`;
    const cmd = [`runme ${S.runopts} -e yelmox`, `-n ${base().file}`, `-o ${outdir}`, ...(out.length ? ["-p " + out.join(" \\\n      ")] : [])].join(" \\\n   ");
    const ta = h("textarea", { class: "fb-ta", rows: 3, placeholder: `paste a -p list, e.g. ctrl.tstep_method=rel forcing.${COEF_LIST_KEY}==alpha,dTo,ramp forcing.terms==lgm,ocean,ramp` });
    body.push(note("Only what differs from the base (or the library group, or the schema default); a library group is used by name. Keys outside the forcing groups (from a preset or a pasted list) pass through."),
      h("div", { class: "fb-fields" },
        h("div", { class: "fb-field" }, h("label", {}, h("code", {}, "runme options")), h("span", { class: "fb-val" }, h("input", { type: "text", class: "fb-text", value: S.runopts, onchange: (e) => { S.runopts = e.target.value; renderRight(); } }))),
        h("div", { class: "fb-field" }, h("label", {}, h("code", {}, "-o")), h("span", { class: "fb-val" }, h("input", { type: "text", class: "fb-text", value: outdir, onchange: (e) => { S.outdir = e.target.value; renderRight(); } })))),
      h("pre", { class: "fb-cfg" }, h("code", {}, cmd)),
      h("div", {}, h("button", { class: "fb-btn", onclick: () => navigator.clipboard.writeText(cmd.replace(/ \\\n\s+/g, " ")) }, "copy (one line)")),
      ...problems.map(warn),
      ...(blocks.length ? [h("pre", { class: "fb-cfg" }, h("code", {}, blocks.join("\n\n")))] : []),
      h("h4", {}, "Load a -p list onto this base"), ta,
      h("div", {}, h("button", { class: "fb-btn", onclick: () => { const errs = parseP(ta.value); S.preset = null; render(); if (errs.length) alert(errs.join("\n")); } }, "load -p list")));
  }
  return h("div", { class: "fb-eqblock fb-output" }, h("h3", {}, "The run, written out"), head, ...body);
}

// --------------------------------------------------------------------- render

// The right column: block 0 (the time slider) and BELOW, the rest; a move of
// the slider redraws BELOW only.
let LEFT = null, RIGHT = null, BELOW = null;
function renderBelow() {
  const tl = timeline();
  const y = RIGHT.scrollTop;
  TAUVAL.textContent = tauText(tl, evalTime(tl));
  BELOW.replaceChildren(...rightEquations(tl), rightRun(tl), rightOutput());
  RIGHT.scrollTop = y;
}
function renderRight() {
  const tl = timeline();
  const y = RIGHT.scrollTop;
  BELOW = BELOW || h("div", { class: "fb-below" });
  RIGHT.replaceChildren(rightTime(tl), BELOW);
  renderBelow();
  RIGHT.scrollTop = y;
}
function render() {
  const y = window.scrollY;
  const tl = timeline();
  root.replaceChildren();
  LEFT = h("div", { class: "fb-left" }, secTime(tl), secEquation(tl), secAfter(),
    h("p", { class: "fb-help" }, `Data: forcing.json generated ${D.generated}${D.commit ? " from " + D.commit : ""} (scripts/forcing-builder/export.jl).`));
  RIGHT = RIGHT || h("div", { class: "fb-right" });
  root.append(explainer(), topBar(), h("div", { class: "fb-cols" }, LEFT, RIGHT));
  renderRight();
  window.scrollTo(0, y);
}

fetch(new URL("forcing-builder/forcing.json", document.baseURI)).then(r => r.json()).then(d => {
  D = d;
  TEMPLATES = (D.templates || []).map(t => [new RegExp(t.pattern), t.template]);
  const q = new URLSearchParams(location.search);
  const pre = (D.presets || []).find(p => p.name === q.get("preset"));
  if (pre) { loadPreset(pre); return; }
  loadBase(D.bases[q.get("base")] ? q.get("base") : (D.bases.yelmox_Greenland ? "yelmox_Greenland" : Object.keys(D.bases)[0]));
  if (q.get("p")) parseP(q.get("p"));
  render();
}).catch(err => { root.replaceChildren(h("p", { class: "fb-warn" }, "forcing.json could not be loaded: " + err)); });
