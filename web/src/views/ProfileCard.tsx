import { forwardRef, memo, useEffect, useImperativeHandle, useMemo, useRef, useState } from "react";
import { fetchCell, type CellDetail, type Meta } from "../api/client";
import { FONT, axes, chrLabel, css, esc, fit, fmt, fmtBp, fmtCount, fmtInt, niceTicks, tip } from "../lib/canvas";
import { barEnd, binAt, binRange, chromLabels, typicalSpan } from "../lib/genome";
import { fails } from "../lib/filters";
import { fmtMetric, isOutlier } from "../lib/metrics";
import { DIVERGING, PALETTES } from "../lib/palettes";
import { discreteColor, isCall, scaleFor } from "../lib/scales";
import { useDerived } from "../state/derived";
import { ALLELE, useStore } from "../state/store";
import { FigureDialog } from "./FigureDialog";
import { DlIcon, useSize, useWidth } from "./widgets";

const PR = { l: 50, t: 12, r: 10, b: 26 };
type Vals = (number | null)[];

// recent cells, shared by every card that shows the selected cell
const cellCache = new WeakMap<Meta, Map<number, Promise<CellDetail>>>();
function cachedCell(meta: Meta, index: number): Promise<CellDetail> {
  let m = cellCache.get(meta);
  if (!m) cellCache.set(meta, (m = new Map()));
  let p = m.get(index);
  if (!p) {
    p = fetchCell(index);
    p.catch(() => m.delete(index));
    m.set(index, p);
    if (m.size > 32) m.delete(m.keys().next().value!);
  }
  return p;
}

/** The selected cell's layers, fetched when the selection changes. */
export function useCellDetail(index: number | null): CellDetail | null {
  const { meta } = useDerived();
  const [cell, setCell] = useState<{ i: number; c: CellDetail } | null>(null);
  useEffect(() => {
    if (index == null) return;
    let live = true;
    cachedCell(meta, index).then((c) => live && setCell({ i: index, c })).catch(() => live && setCell(null));
    return () => { live = false; };
  }, [index, meta]);
  return index != null && cell?.i === index ? cell.c : null;
}

function drawChromosomes(g: CanvasRenderingContext2D, meta: Meta, xS: (v: number) => number, gx0: number, gx1: number, top: number, h: number, labelY: number | null) {
  g.textAlign = "center"; g.textBaseline = "alphabetic"; g.font = `11px ${FONT}`;
  const spans: [string, number, number][] = [];
  for (const { chr, start_cum, end_cum } of meta.chr_info) {
    if (end_cum < gx0 || start_cum > gx1) continue;
    const a = xS(Math.max(start_cum, gx0)), b = xS(Math.min(end_cum, gx1));
    if (start_cum > gx0) { g.strokeStyle = css("--line2"); g.lineWidth = 1; g.beginPath(); g.moveTo(Math.round(a) + 0.5, top); g.lineTo(Math.round(a) + 0.5, top + h); g.stroke(); }
    spans.push([chr, a, b]);
  }
  if (labelY == null) return;
  g.fillStyle = css("--ink2");
  for (const [t, x] of chromLabels(g, spans)) g.fillText(t, x, labelY);
}

/** Runs of one state within a chromosome: [first bin, last bin + 1, state]. */
function runsOf(chr: string[], vals: Vals): [number, number, number][] {
  const out: [number, number, number][] = [], nb = vals.length;
  for (let s = 0; s < nb; ) {
    let e = s + 1;
    while (e < nb && vals[e] === vals[s] && chr[e] === chr[s]) e++;
    if (vals[s] != null) out.push([s, e, vals[s] as number]);
    s = e;
  }
  return out;
}

/* ── profile: copy-number calls over the bins, following the heatmap's genomic zoom ── */
const ProfilePlot = memo(forwardRef<HTMLCanvasElement, { cell: CellDetail; height?: number }>(function ProfilePlot({ cell, height = 300 }, ref) {
  const d = useDerived();
  const { meta, totalLayer, binsLayer, sem } = d;
  const view = useStore((s) => s.view);
  const zoom = useStore((s) => s.zoom);
  const showCi = useStore((s) => s.showCi);
  const theme = useStore((s) => s.theme);
  const box = useRef<HTMLDivElement>(null);
  const cv = useRef<HTMLCanvasElement>(null);
  useImperativeHandle(ref, () => cv.current!, []);
  const width = useWidth(box);
  const B = meta.bins, nb = meta.n_bins, genomeMax = B.end_cum[nb - 1];
  const [gx0, gx1] = zoom?.x ?? [0, genomeMax];

  const L = cell.layers;
  const isAS = view === ALLELE && !!L.major_cn && !!L.minor_cn;
  // calls: the integer total, or the view's own call layer when there is none
  const callLayer = totalLayer ?? (isCall(sem[view] ?? "count") ? view : null);
  const prof = callLayer ? L[callLayer] : null;
  const raw = useMemo(() => {
    if (!binsLayer || !L[binsLayer[0]]) return !callLayer && L[view] ? L[view] : null;
    const v = L[binsLayer[0]];
    if (binsLayer[1] === "cn" || !prof) return v;
    let s = 0, n = 0;
    for (const x of prof) if (x != null) { s += x; n++; }
    const psi = n ? s / n : 2;
    return v.map((x) => (x == null ? null : psi * Math.pow(2, x)));
  }, [L, binsLayer, prof, callLayer, view]);
  const ivMaj = meta.intervals.major_cn, ivMin = meta.intervals.minor_cn;
  const hasCi = isAS && showCi && !!ivMaj && !!ivMin && L[ivMaj[0]]?.some((v) => v != null);
  const cnScale = !!callLayer && (sem[callLayer] === "integer_cn" || sem[callLayer] === "allele_cn");
  const ternary = !!callLayer && sem[callLayer] === "ternary_state";
  const palette = useStore((s) => s.palette);
  const stateScale = useMemo(() => scaleFor("integer_cn", { palette: PALETTES[palette], dark: theme === "dark" }), [palette, theme]);
  const ternaryScale = useMemo(() => scaleFor("ternary_state", { palette: PALETTES[palette], dark: theme === "dark" }), [palette, theme]);

  useEffect(() => {
    const c = cv.current;
    if (!c) return;
    const g = fit(c, width, height);
    g.clearRect(0, 0, width, height);
    const plotW = width - PR.l - PR.r, plotH = height - PR.t - PR.b;
    const xS = (v: number) => PR.l + ((v - gx0) / (gx1 - gx0)) * plotW;
    const [b0, b1] = binRange(B.start_cum, B.end_cum, gx0, gx1);
    const nz = b1 - b0;

    let Y: (v: number) => number, ymax = 0;
    if (cnScale) {
      const top = (arr: Vals) => { let m = 0; for (const v of arr) if (v != null && v > m) m = v; return m; };
      ymax = Math.max(5, Math.min(12, (isAS ? top(L.major_cn) : top(prof!)) + 1));
      Y = (v) => PR.t + plotH - (Math.max(-0.3, Math.min(ymax + 0.4, v)) / (ymax + 0.4)) * plotH;
      g.font = `11px ${FONT}`; g.strokeStyle = css("--grid"); g.lineWidth = 1; g.fillStyle = css("--ink3"); g.textAlign = "right"; g.textBaseline = "middle";
      for (let y = 0; y <= ymax; y++) {
        g.setLineDash(y === 2 ? [] : [3, 4]); g.beginPath(); g.moveTo(PR.l, Y(y)); g.lineTo(PR.l + plotW, Y(y)); g.stroke();
        g.fillText(String(y), PR.l - 8, Y(y));
      }
      g.setLineDash([]);
    } else if (ternary) {
      // gain / normal / loss calls (epiAneufinder)
      Y = (v) => PR.t + plotH - ((Math.max(-1.5, Math.min(1.5, v)) + 1.5) / 3) * plotH;
      g.font = `11px ${FONT}`; g.strokeStyle = css("--grid"); g.lineWidth = 1; g.fillStyle = css("--ink3"); g.textAlign = "right"; g.textBaseline = "middle";
      ["Loss", "Normal", "Gain"].forEach((t, k) => {
        g.setLineDash(k === 1 ? [] : [3, 4]); g.beginPath(); g.moveTo(PR.l, Y(k - 1)); g.lineTo(PR.l + plotW, Y(k - 1)); g.stroke();
        g.fillText(t, PR.l - 8, Y(k - 1));
      });
      g.setLineDash([]);
    } else {
      // a generic layer: its own range
      const vals = (raw ?? prof ?? []).slice(b0, b1).filter((v): v is number => v != null).sort((a, b) => a - b);
      let lo = vals.length ? vals[Math.floor(0.005 * (vals.length - 1))] : 0, hi = vals.length ? vals[Math.floor(0.995 * (vals.length - 1))] : 1;
      if (prof) for (let i = b0; i < b1; i++) { const v = prof[i]; if (v != null) { lo = Math.min(lo, v); hi = Math.max(hi, v); } }
      if (!(hi > lo)) { lo -= 1; hi += 1; }
      const pad = (hi - lo) * 0.06; lo -= pad; hi += pad;
      Y = (v) => PR.t + plotH - ((Math.max(lo, Math.min(hi, v)) - lo) / (hi - lo)) * plotH;
      g.font = `11px ${FONT}`; g.strokeStyle = css("--grid"); g.lineWidth = 1; g.fillStyle = css("--ink3"); g.textAlign = "right"; g.textBaseline = "middle";
      for (const t of niceTicks(lo, hi, 5)) { g.beginPath(); g.moveTo(PR.l, Y(t)); g.lineTo(PR.l + plotW, Y(t)); g.stroke(); g.fillText(String(t), PR.l - 8, Y(t)); }
    }
    drawChromosomes(g, meta, xS, gx0, gx1, PR.t, plotH, PR.t + plotH + 18);
    g.save(); g.beginPath(); g.rect(PR.l, 0, plotW, height); g.clip();
    if (raw) {
      // bins coloured by their copy-number state, as in scgenome / SIGNALS profiles
      const pointCol = css("--points");
      g.globalAlpha = 0.8;
      const r = nz > 1500 ? 1.4 : nz > 400 ? 1.9 : 2.7;
      for (let i = b0; i < b1; i++) {
        const v = raw[i];
        if (v == null || !Number.isFinite(v)) continue;
        g.fillStyle = (cnScale && prof && discreteColor(stateScale, prof[i] ?? NaN)) || pointCol;
        g.beginPath(); g.arc(xS((B.start_cum[i] + B.end_cum[i]) / 2), Y(v), r, 0, 7); g.fill();
      }
      g.globalAlpha = 1;
    }
    const cap = (v: number) => (cnScale ? Math.min(v, ymax) : v);
    const runs = prof ? runsOf(B.chr, prof) : [];
    if (prof && ternary) {
      g.lineWidth = 4; g.lineCap = "round";
      for (const [s, e, v] of runs) {
        g.strokeStyle = discreteColor(ternaryScale, v) || css("--call");
        g.beginPath(); g.moveTo(xS(B.start_cum[s]) + 1, Y(v)); g.lineTo(xS(B.end_cum[e - 1]) - 1, Y(v)); g.stroke();
      }
    } else if (prof && !isAS) {
      g.strokeStyle = css("--call"); g.lineWidth = 3; g.lineCap = "round"; g.beginPath();
      for (const [s, e, v] of runs) { g.moveTo(xS(B.start_cum[s]) + 1, Y(cap(v))); g.lineTo(xS(B.end_cum[e - 1]) - 1, Y(cap(v))); }
      g.stroke();
    } else if (isAS) {
      const sets: [Vals, string, [string, string] | undefined, number][] = [[L.major_cn, "--major", ivMaj, 0.07], [L.minor_cn, "--minor", ivMin, -0.07]];
      for (const [arr, col, iv, off] of sets) {
        if (hasCi && iv) {
          const lo = L[iv[0]], hi = L[iv[1]];
          g.fillStyle = css(col); g.globalAlpha = 0.2;
          for (let i = b0; i < b1; i++) {
            const a = arr[i], l = lo[i], h = hi[i];
            if (a == null || l == null || h == null) continue;
            const yt = Y(Math.max(a, h)), yb = Y(Math.min(a, l));
            g.fillRect(xS(B.start_cum[i]), yt, xS(B.end_cum[i]) - xS(B.start_cum[i]), yb - yt);
          }
          g.globalAlpha = 1;
        }
        g.strokeStyle = css(col); g.lineWidth = 3; g.lineCap = "round"; g.beginPath();
        for (const [s, e, v] of runsOf(B.chr, arr)) { g.moveTo(xS(B.start_cum[s]) + 1, Y(Math.min(v, ymax) + off)); g.lineTo(xS(B.end_cum[e - 1]) - 1, Y(Math.min(v, ymax) + off)); }
        g.stroke();
      }
    }
    if (raw && prof) {
      g.strokeStyle = css("--mean"); g.lineWidth = 1.2; g.lineCap = "butt"; g.beginPath();
      for (const [s, e] of runs) {
        let m = 0, k = 0;
        for (let i = s; i < e; i++) { const v = raw[i]; if (v != null && Number.isFinite(v)) { m += v; k++; } }
        if (!k) continue;
        m /= k;
        g.moveTo(xS(B.start_cum[s]), Y(m)); g.lineTo(xS(B.end_cum[e - 1]), Y(m));
      }
      g.stroke();
    }
    g.restore();
    if (!raw && prof && !ternary && !d.exprLayer) { g.fillStyle = css("--ink3"); g.textAlign = "left"; g.textBaseline = "alphabetic"; g.font = `11px ${FONT}`; g.fillText("no per-bin values in this dataset; showing the calls", PR.l + 8, PR.t + 14); }
  }, [cell, width, height, gx0, gx1, isAS, raw, prof, hasCi, cnScale, L, meta, B, nb, ivMaj, ivMin, theme, stateScale, ternary, ternaryScale, d.exprLayer]);

  const move = (e: React.MouseEvent) => {
    const r = cv.current!.getBoundingClientRect(), x = e.clientX - r.left, plotW = width - PR.l - PR.r;
    if (x < PR.l || x > PR.l + plotW) return tip(null);
    const b = binAt(B.end_cum, gx0 + ((x - PR.l) / plotW) * (gx1 - gx0));
    const line = isAS
      ? `major|minor <b>${L.major_cn[b] ?? "–"}|${L.minor_cn[b] ?? "–"}</b>`
      : prof ? `${esc(callLayer === totalLayer ? "copy number" : callLayer!)} <b>${prof[b] ?? "–"}</b>` : "";
    const rd = d.countLayer ? L[d.countLayer]?.[b] : null;
    const hp = isAS && d.bafLayer && d.bafNLayer && L[d.bafLayer] && L[d.bafNLayer] ? hapReads(L[d.bafLayer], L[d.bafNLayer], b) : null;
    tip(e, `<span class="k">${chrLabel(B.chr[b])}:${fmtBp(B.start[b])}–${fmtBp(B.end[b])}</span>${line ? "<br>" + line : ""}` +
      (raw && raw[b] != null && !cnScale ? `<br><span class="k">${esc(view)}</span> ${fmt(raw[b])}` : "") +
      (rd != null ? `<br><span class="k">reads</span> ${rd}` : "") +
      (hp ? `<br><span class="k">hap1 | hap2 reads</span> ${hp[0]} | ${hp[1]}` : ""));
  };

  const items: [string, string, string][] = [];
  if (raw) items.push(["dot", cnScale ? "--points-key" : "--points", cnScale ? "bin, on the copy-number scale, coloured by its state" : `bin (${view})`]);
  if (raw && prof) items.push(["line", "--mean", "segment mean"]);
  if (isAS) items.push(["thick", "--major", "major allele"], ["thick", "--minor", "minor allele"]);
  else if (prof) items.push(["thick", "--call", callLayer === totalLayer ? "integer copy number" : callLayer!]);
  return (
    <div ref={box} style={{ width: "100%" }}>
      <canvas ref={cv} onMouseMove={move} onMouseLeave={() => tip(null)} />
      <div className="legend lines">
        {items.map(([k, c, t]) => (
          <span className="sw" key={t}><i style={{ width: k === "dot" ? 8 : 20, height: k === "thick" ? 4 : k === "dot" ? 8 : 2, borderRadius: k === "dot" ? 4 : 2, background: `var(${c})` }} />{t}</span>
        ))}
        {hasCi && <span className="sw"><i style={{ width: 20, height: 8, borderRadius: 2, background: "var(--major)", opacity: 0.25 }} />BAF 5–95% interval</span>}
      </div>
    </div>
  );
}));

/* ── BAF: per-bin allele fraction, with the BAF the allele-specific calls expect ── */
function BafPlot({ baf, n, major, minor, height = 180 }: { baf: Vals; n: Vals | null; major: Vals | null; minor: Vals | null; height?: number }) {
  const { meta } = useDerived();
  const zoom = useStore((s) => s.zoom);
  const theme = useStore((s) => s.theme);
  const box = useRef<HTMLDivElement>(null);
  const cv = useRef<HTMLCanvasElement>(null);
  const width = useWidth(box);
  const B = meta.bins, genomeMax = B.end_cum[meta.n_bins - 1];
  const [gx0, gx1] = zoom?.x ?? [0, genomeMax];
  const T = 18, BOT = 6;
  useEffect(() => {
    const c = cv.current;
    if (!c) return;
    const g = fit(c, width, height);
    g.clearRect(0, 0, width, height);
    const plotW = width - PR.l - PR.r, plotH = height - T - BOT;
    const xS = (v: number) => PR.l + ((v - gx0) / (gx1 - gx0)) * plotW;
    const Y = (v: number) => T + plotH - Math.max(0, Math.min(1, v)) * plotH;
    const [b0, b1] = binRange(B.start_cum, B.end_cum, gx0, gx1);
    g.font = `11px ${FONT}`; g.lineWidth = 1; g.fillStyle = css("--ink3"); g.textAlign = "right"; g.textBaseline = "middle";
    for (const [v, t] of [[0, "0"], [0.2, "0.2"], [0.4, "0.4"], [0.5, "0.5"], [0.6, "0.6"], [0.8, "0.8"], [1, "1"]] as const) {
      g.strokeStyle = css(v === 0.5 ? "--ink3" : "--grid"); g.setLineDash(v === 0.5 ? [5, 4] : []);
      g.beginPath(); g.moveTo(PR.l, Y(v)); g.lineTo(PR.l + plotW, Y(v)); g.stroke();
      g.fillText(t, PR.l - 8, Y(v));
    }
    g.setLineDash([]);
    drawChromosomes(g, meta, xS, gx0, gx1, T, plotH, null);
    g.save(); g.beginPath(); g.rect(PR.l, 0, plotW, height); g.clip();
    // point size from the evidence (reads or SNPs), against the cell's own 95th percentile
    const ev: number[] = [];
    if (n) for (let i = b0; i < b1; i++) { const v = n[i]; if (v != null && v > 0) ev.push(v); }
    ev.sort((a, b) => a - b);
    const top = ev.length ? ev[Math.floor(0.95 * (ev.length - 1))] : 1;
    g.fillStyle = css("--series1"); g.globalAlpha = 0.6;
    for (let i = b0; i < b1; i++) {
      const v = baf[i];
      if (v == null || !Number.isFinite(v)) continue;
      const w = n ? n[i] ?? 0 : 1;
      const r = n ? 1.2 + 2.8 * Math.sqrt(Math.min(1, w / Math.max(1, top))) : 1.8;
      g.beginPath(); g.arc(xS((B.start_cum[i] + B.end_cum[i]) / 2), Y(v), r, 0, 7); g.fill();
    }
    g.globalAlpha = 1;
    // dotted: the BAF the allele-specific state predicts, minor / total and major / total
    if (major && minor) {
      const key: Vals = major.map((a, i) => (a == null || minor[i] == null ? null : (a as number) * 1000 + (minor[i] as number)));
      g.strokeStyle = css("--ink"); g.lineWidth = 1.4; g.setLineDash([2, 3]); g.beginPath();
      for (const [s, e, k] of runsOf(B.chr, key)) {
        const a = Math.floor(k / 1000), b = k % 1000, t = a + b;
        if (!(t > 0) || a === b) continue;
        for (const f of [b / t, a / t]) { g.moveTo(xS(B.start_cum[s]), Y(f)); g.lineTo(xS(B.end_cum[e - 1]), Y(f)); }
      }
      g.stroke(); g.setLineDash([]);
    }
    g.restore();
    g.fillStyle = css("--ink2"); g.textAlign = "left"; g.textBaseline = "alphabetic";
    g.fillText("BAF per bin", PR.l, 12);
    g.fillStyle = css("--ink3"); g.textAlign = "right";
    g.fillText(major && minor ? "dotted: BAF expected from the allele-specific calls" + (n ? " · size: reads" : "") : n ? "size: reads" : "", PR.l + plotW, 12);
  }, [baf, n, major, minor, width, height, gx0, gx1, meta, B, theme]);
  const move = (e: React.MouseEvent<HTMLCanvasElement>) => {
    const b = binUnder(e, width, gx0, gx1, B.end_cum);
    tip(b == null ? null : e, b == null ? "" : bafTip(meta, b, baf, n, major, minor));
  };
  return <div ref={box} style={{ width: "100%", marginTop: 10 }}><canvas ref={cv} onMouseMove={move} onMouseLeave={() => tip(null)} /></div>;
}

/** hap1 / hap2 reads of a bin from its BAF (hap2 / (hap1 + hap2)) and its reads, or null. */
function hapReads(baf: Vals, n: Vals | null, i: number): [number, number] | null {
  const b = baf[i], t = n ? n[i] : null;
  if (b == null || t == null || !Number.isFinite(b) || !(t > 0)) return null;
  const h2 = Math.round(b * t);
  return [t - h2, h2];
}

/** Tooltip for a bin of the BAF panel: position, BAF, hap1 | hap2 reads
 *  and the call with the BAF it predicts. */
function bafTip(meta: Meta, b: number, baf: Vals, n: Vals | null, major: Vals | null, minor: Vals | null): string {
  const B = meta.bins;
  const v = baf[b], h = hapReads(baf, n, b);
  let out = `<span class="k">${chrLabel(B.chr[b])}:${fmtBp(B.start[b])}–${fmtBp(B.end[b])}</span>`;
  out += `<br><span class="k">BAF</span> <b>${v != null && Number.isFinite(v) ? fmt(v, 2) : "–"}</b>`;
  if (h) out += `<br><span class="k">hap1 | hap2 reads</span> <b>${h[0]} | ${h[1]}</b>`;
  else if (n && n[b] != null) out += `<br><span class="k">reads</span> ${n[b]}`;
  const a = major?.[b], c = minor?.[b];
  if (a != null && c != null) {
    const t = a + c;
    out += `<br><span class="k">major|minor</span> ${a}|${c}` + (t > 0 ? ` <span class="k">· expected BAF</span> ${a === c ? "0.5" : `${fmt(c / t, 2)} or ${fmt(a / t, 2)}`}` : "");
  }
  return out;
}

/** The bin under the mouse in a panel drawn with the profile's margins and zoom. */
function binUnder(e: React.MouseEvent<HTMLCanvasElement>, width: number, gx0: number, gx1: number, endCum: ArrayLike<number>): number | null {
  const r = e.currentTarget.getBoundingClientRect(), x = e.clientX - r.left, plotW = width - PR.l - PR.r;
  if (x < PR.l || x > PR.l + plotW) return null;
  return binAt(endCum, gx0 + ((x - PR.l) / plotW) * (gx1 - gx0));
}

/* ── expression (scRNA / multiome): relative expression per bin, on its own axis ── */
/** The gain / normal / loss layer of a cell, if the dataset has one. */
function ternaryOf(L: Record<string, Vals>, sem: Record<string, string>): Vals | null {
  const k = Object.keys(sem).find((n) => sem[n] === "ternary_state" && L[n]);
  return k ? L[k] : null;
}

/** "its own cell_type (16 with >= 20 cells; the others: all cells)" -> "its own cell_type" */
const refShort = (r: string) => r.replace(/\s*\(.*\)$/, "");

function ExpressionPlot({ expr, call, callName, height = 170 }: { expr: Vals; call: Vals | null; callName: string; height?: number }) {
  const { meta } = useDerived();
  // scATAC: logR of fragment counts against the same bin in reference cells (epiAneufinder)
  const atac = meta.modality === "scATAC";
  const zoom = useStore((s) => s.zoom);
  const theme = useStore((s) => s.theme);
  const box = useRef<HTMLDivElement>(null);
  const cv = useRef<HTMLCanvasElement>(null);
  const width = useWidth(box);
  const B = meta.bins, genomeMax = B.end_cum[meta.n_bins - 1];
  const [gx0, gx1] = zoom?.x ?? [0, genomeMax];
  const T = 18, BOT = 6;
  // the axis: ± the 98th percentile of this cell's |log2 expression|, at least ±0.05
  const range = useMemo(() => {
    const v: number[] = [];
    for (const x of expr) if (x != null && Number.isFinite(x)) v.push(Math.abs(x));
    v.sort((a, b) => a - b);
    return v.length ? Math.max(0.05, v[Math.floor(0.98 * (v.length - 1))] * 1.1) : 0.2;
  }, [expr]);
  useEffect(() => {
    const c = cv.current;
    if (!c) return;
    const g = fit(c, width, height);
    g.clearRect(0, 0, width, height);
    const plotW = width - PR.l - PR.r, plotH = height - T - BOT;
    const xS = (v: number) => PR.l + ((v - gx0) / (gx1 - gx0)) * plotW;
    const Y = (v: number) => T + plotH / 2 - (Math.max(-range, Math.min(range, v)) / range) * (plotH / 2);
    const [b0, b1] = binRange(B.start_cum, B.end_cum, gx0, gx1);
    // ticks in relative expression (×), placed on the log2 axis; scATAC: logR itself
    g.font = `11px ${FONT}`; g.lineWidth = 1; g.fillStyle = css("--ink3"); g.textAlign = "right"; g.textBaseline = "middle";
    const ticks: [number, string][] = [];
    if (atac) {
      const st = range > 2.5 ? 1 : 0.5;
      for (let v = -Math.floor(range / st) * st; v <= range + 1e-9; v += st) ticks.push([v, v === 0 ? "0" : (v > 0 ? "+" : "") + v]);
    } else {
      const step = range > 0.4 ? 0.25 : range > 0.15 ? 0.1 : 0.05;
      ticks.push([0, "1×"]);
      for (let k = 1; k < 20; k++) {
        const up = 1 + k * step, dn = 1 - k * step;
        if (Math.log2(up) <= range) ticks.push([Math.log2(up), `${Math.round(up * 100) / 100}×`]);
        if (dn > 0 && -Math.log2(dn) <= range) ticks.push([Math.log2(dn), `${Math.round(dn * 100) / 100}×`]);
      }
    }
    for (const [v, label] of ticks) {
      const y = Y(v);
      g.strokeStyle = css(v === 0 ? "--ink3" : "--grid"); g.setLineDash(v === 0 ? [5, 4] : []);
      g.beginPath(); g.moveTo(PR.l, y); g.lineTo(PR.l + plotW, y); g.stroke();
      g.fillText(label, PR.l - 6, y);
    }
    g.setLineDash([]);
    drawChromosomes(g, meta, xS, gx0, gx1, T, plotH, null);
    g.save(); g.beginPath(); g.rect(PR.l, 0, plotW, height); g.clip();
    const lo = css("--expr-lo") || DIVERGING[0], hi = css("--expr-hi") || DIVERGING[4], mid = css("--ink3");
    for (let i = b0; i < b1; i++) {
      const v = expr[i];
      if (v == null || !Number.isFinite(v)) continue;
      const f = Math.min(1, Math.abs(v) / range);
      g.fillStyle = f < 0.15 ? mid : v > 0 ? hi : lo;
      g.globalAlpha = 0.35 + 0.6 * f;
      g.beginPath(); g.arc(xS((B.start_cum[i] + B.end_cum[i]) / 2), Y(v), 1.8, 0, 7); g.fill();
    }
    g.globalAlpha = 1;
    if (atac && call) {
      // single bins are noisy: the mean logR over each run of one call
      g.strokeStyle = css("--mean"); g.lineWidth = 1.4; g.beginPath();
      for (const [s, e] of runsOf(B.chr, call)) {
        if (e <= b0 || s >= b1) continue;
        let m = 0, k = 0;
        for (let i = s; i < e; i++) { const v = expr[i]; if (v != null && Number.isFinite(v)) { m += v; k++; } }
        if (!k) continue;
        m /= k;
        g.moveTo(xS(B.start_cum[s]), Y(m)); g.lineTo(xS(B.end_cum[e - 1]), Y(m));
      }
      g.stroke();
    }
    g.restore();
    g.textBaseline = "alphabetic"; g.textAlign = "left"; g.fillStyle = css("--ink2");
    g.fillText(atac ? `logR per bin (fragments, against the same bin in ${meta.logr_reference ? refShort(meta.logr_reference) : "all cells"})` : "Expression per bin (relative to the reference cells)", PR.l, 12);
    g.textAlign = "right"; g.fillStyle = css("--ink3");
    g.fillText(atac ? "line: mean over each gain / normal / loss stretch" : "not on the copy-number scale: expression is compressed and smoothed", PR.l + plotW, 12);
  }, [expr, call, atac, range, width, height, gx0, gx1, meta, B, theme]);
  const move = (e: React.MouseEvent<HTMLCanvasElement>) => {
    const b = binUnder(e, width, gx0, gx1, B.end_cum);
    if (b == null) return tip(null);
    const v = expr[b];
    tip(e, `<span class="k">${chrLabel(B.chr[b])}:${fmtBp(B.start[b])}–${fmtBp(B.end[b])}</span>` +
      (atac
        ? `<br><span class="k">logR</span> <b>${v != null && Number.isFinite(v) ? fmt(v, 2) : "–"}</b>`
        : `<br><span class="k">expression</span> <b>${v != null && Number.isFinite(v) ? fmt(2 ** v, 3) + "×" : "–"}</b>`) +
      (call && call[b] != null ? `<br><span class="k">${esc(callName)}</span> ${atac ? (["loss", "normal", "gain"][call[b]! + 1] ?? call[b]) : call[b]}` : ""));
  };
  return <div ref={box} style={{ width: "100%", marginTop: 10 }}><canvas ref={cv} onMouseMove={move} onMouseLeave={() => tip(null)} /></div>;
}

/* ── read depth: per-bin read counts on the profile's x scale ── */
function ReadDepthPlot({ reads, height = 120 }: { reads: Vals; height?: number }) {
  const { meta } = useDerived();
  const zoom = useStore((s) => s.zoom);
  const theme = useStore((s) => s.theme);
  const box = useRef<HTMLDivElement>(null);
  const cv = useRef<HTMLCanvasElement>(null);
  const width = useWidth(box);
  const B = meta.bins, genomeMax = B.end_cum[meta.n_bins - 1];
  const [gx0, gx1] = zoom?.x ?? [0, genomeMax];
  const T = 18, BOT = 6;
  const typical = useMemo(() => typicalSpan(B.start_cum, B.end_cum), [B]);
  useEffect(() => {
    const c = cv.current;
    if (!c) return;
    const g = fit(c, width, height);
    g.clearRect(0, 0, width, height);
    const plotW = width - PR.l - PR.r, plotH = height - T - BOT;
    const xS = (v: number) => PR.l + ((v - gx0) / (gx1 - gx0)) * plotW;
    const [b0, b1] = binRange(B.start_cum, B.end_cum, gx0, gx1);
    // the axis stops at the 99.5th percentile so a few extreme bins don't flatten the rest;
    // clipped bars get a mark on top
    let sum = 0, n = 0;
    const seen: number[] = [];
    for (let i = b0; i < b1; i++) { const v = reads[i]; if (v == null) continue; seen.push(v); sum += v; n++; }
    seen.sort((a, b) => a - b);
    const max = Math.max(1, seen.length ? seen[Math.floor(0.995 * (seen.length - 1))] * 1.1 : 1);
    const Y = (v: number) => T + plotH - (Math.min(v, max) / max) * plotH;
    g.font = `11px ${FONT}`; g.strokeStyle = css("--grid"); g.lineWidth = 1; g.fillStyle = css("--ink3"); g.textAlign = "right"; g.textBaseline = "middle";
    for (const t of niceTicks(0, max, 3)) { g.beginPath(); g.moveTo(PR.l, Y(t)); g.lineTo(PR.l + plotW, Y(t)); g.stroke(); g.fillText(String(t), PR.l - 8, Y(t)); }
    drawChromosomes(g, meta, xS, gx0, gx1, T, plotH, null);
    g.save(); g.beginPath(); g.rect(PR.l, 0, plotW, height); g.clip();
    // a bin much wider than typical (a centromere or other unmappable stretch merged into it):
    // its bar is a typical bin wide, and the rest of its span is hatched, with no reads drawn
    g.strokeStyle = css("--grid"); g.lineWidth = 1;
    for (let i = b0; i < b1; i++) {
      const e = barEnd(B.start_cum, B.end_cum, i, typical);
      if (e >= B.end_cum[i]) continue;
      const a = xS(e), b = xS(B.end_cum[i]);
      g.save(); g.beginPath(); g.rect(a, T, b - a, plotH); g.clip();
      g.beginPath();
      for (let x = a - plotH; x < b; x += 6) { g.moveTo(x, T + plotH); g.lineTo(x + plotH, T); }
      g.stroke(); g.restore();
    }
    g.fillStyle = css("--series1"); g.globalAlpha = 0.75;
    for (let i = b0; i < b1; i++) {
      const v = reads[i];
      if (v == null) continue;
      const x0 = xS(B.start_cum[i]), x1 = xS(barEnd(B.start_cum, B.end_cum, i, typical)), y = Y(v);
      g.fillRect(x0, y, Math.max(0.6, x1 - x0), T + plotH - y);
    }
    g.globalAlpha = 1;
    g.fillStyle = css("--serious");
    for (let i = b0; i < b1; i++) {
      const v = reads[i];
      if (v != null && v > max) g.fillRect(xS(B.start_cum[i]), T - 2, Math.max(1.5, xS(barEnd(B.start_cum, B.end_cum, i, typical)) - xS(B.start_cum[i])), 2);
    }
    if (n) { g.strokeStyle = css("--ink3"); g.setLineDash([4, 4]); g.beginPath(); g.moveTo(PR.l, Y(sum / n)); g.lineTo(PR.l + plotW, Y(sum / n)); g.stroke(); g.setLineDash([]); }
    g.restore();
    g.fillStyle = css("--ink2"); g.textAlign = "left"; g.textBaseline = "alphabetic";
    g.fillText("Reads per bin", PR.l, 12);
    if (n) { g.fillStyle = css("--ink3"); g.textAlign = "right"; g.fillText(`dashed: mean ${fmt(sum / n, 1)}`, PR.l + plotW, 12); }
  }, [reads, width, height, gx0, gx1, meta, B, theme, typical]);
  const move = (e: React.MouseEvent<HTMLCanvasElement>) => {
    const b = binUnder(e, width, gx0, gx1, B.end_cum);
    if (b == null) return tip(null);
    const v = reads[b];
    const span = B.end[b] - B.start[b];
    const wide = barEnd(B.start_cum, B.end_cum, b, typical) < B.end_cum[b];
    tip(e, `<span class="k">${chrLabel(B.chr[b])}:${fmtBp(B.start[b])}–${fmtBp(B.end[b])}</span>` +
      `<br><span class="k">reads</span> <b>${v != null ? fmtInt(v) : "–"}</b>` +
      (wide
        ? `<br><span class="k">one bin over ${fmtBp(span)}: a centromere or other stretch with little<br>mappable sequence was merged into it. Its bar is drawn a typical<br>bin wide; the hatched rest has no reads of its own.</span>`
        : ""));
  };
  return <div ref={box} style={{ width: "100%", marginTop: 10 }}><canvas ref={cv} onMouseMove={move} onMouseLeave={() => tip(null)} /></div>;
}

/* ── drag across any of the cell's panels to zoom them (and the heatmap) to that stretch ── */
function ZoomBrush({ children }: { children: React.ReactNode }) {
  const { meta, shown } = useDerived();
  const zoom = useStore((s) => s.zoom);
  const set = useStore((s) => s.set);
  const box = useRef<HTMLDivElement>(null);
  const width = useWidth(box);
  const [band, setBand] = useState<[number, number] | null>(null);
  const start = useRef<number | null>(null);
  const genomeMax = meta.bins.end_cum[meta.n_bins - 1];
  const [gx0, gx1] = zoom?.x ?? [0, genomeMax];
  const plotW = width - PR.l - PR.r;
  const px = (e: React.MouseEvent) => {
    const x = e.clientX - box.current!.getBoundingClientRect().left;
    return Math.max(PR.l, Math.min(PR.l + plotW, x));
  };
  const down = (e: React.MouseEvent) => {
    if (e.button !== 0) return;
    const x = e.clientX - box.current!.getBoundingClientRect().left;
    if (x < PR.l || x > PR.l + plotW) return;
    start.current = x; setBand([x, x]); e.preventDefault();
  };
  const move = (e: React.MouseEvent) => { if (start.current != null) { setBand([start.current, px(e)]); tip(null); } };
  const up = (e: React.MouseEvent) => {
    if (start.current == null) return;
    const a = Math.min(start.current, px(e)), b = Math.max(start.current, px(e));
    start.current = null; setBand(null);
    if (b - a < 6) return;
    const g = (x: number) => gx0 + ((x - PR.l) / plotW) * (gx1 - gx0);
    set({ zoom: { x: [Math.max(0, g(a)), Math.min(genomeMax, g(b))], y: zoom?.y ?? [0, shown.length] } });
  };
  const reset = () => { if (zoom) set({ zoom: { x: [0, genomeMax], y: zoom.y } }); };
  return (
    <div ref={box} className="brush" onMouseDown={down} onMouseMove={move} onMouseUp={up}
      onMouseLeave={() => { start.current = null; setBand(null); }} onDoubleClick={reset}>
      {children}
      {band && <div className="brushband" style={{ left: Math.min(...band), width: Math.abs(band[1] - band[0]) }} />}
    </div>
  );
}

/* ── the card: status, metric pills, profile and read depth ── */
export function ProfileCard() {
  const d = useDerived();
  const { meta, vals, metAll, passing, reasons, countLayer } = d;
  const selected = useStore((s) => s.selected);
  const view = useStore((s) => s.view);
  const zoom = useStore((s) => s.zoom);
  const set = useStore((s) => s.set);
  const cell = useCellDetail(selected);
  const plot = useRef<HTMLCanvasElement>(null);
  const B = meta.bins, genomeMax = B.end_cum[meta.n_bins - 1];
  const xZoomed = !!zoom && (zoom.x[0] > 0 || zoom.x[1] < genomeMax);
  const zoomText = () => {
    const a = binAt(B.end_cum, zoom!.x[0]), b = binAt(B.end_cum, zoom!.x[1]);
    return `${chrLabel(B.chr[a])}:${fmtBp(B.start[a])} – ${chrLabel(B.chr[b])}:${fmtBp(B.end[b])}`;
  };
  const [fig, setFig] = useState(false);
  const name = selected != null ? meta.cells[selected] : null;
  const valOf = (k: string) => (selected != null && vals[k] ? vals[k][selected] : NaN);
  const failsFilter = (metric: string) => fails(d.filters[metric], valOf(metric));
  const metricOf = (k: string) => metAll.find((m) => m.key === k);

  const pills: [string, string, string, string?][] = [];
  if (selected != null) {
    const add = (k: string, label: string, v: string, title?: string) => { if (!Number.isNaN(valOf(k))) pills.push([k, label, v, title]); };
    add("ploidy", "ploidy", fmt(valOf("ploidy")));
    add("n_reads_total", "reads", fmtCount(valOf("n_reads_total")), `${fmtInt(valOf("n_reads_total"))} reads in this cell's bins`);
    add("reads_per_bin", "reads / bin", fmt(valOf("reads_per_bin"), 1));
    add("coverage", "coverage", fmtMetric("coverage", valOf("coverage")));
    add("mapd", "MAPD", fmt(valOf("mapd"), 3));
    add("median_residual", "residual", fmt(valOf("median_residual"), 3));
    add("gini", "Gini", fmt(valOf("gini"), 3));
    add("bins_with_cna", "CNA bins", fmtInt(valOf("bins_with_cna")), `bins whose copy number is not 2, of ${meta.n_bins}`);
    add("breakpoints", "breakpoints", fmtInt(valOf("breakpoints")));
  }
  const ct = selected != null ? meta.obs.cell_type?.[selected] : null;
  const status = selected == null ? null : passing[selected]
    ? <span className="status" style={{ color: "var(--good)" }}><svg width="12" height="12" viewBox="0 0 12 12"><circle cx="6" cy="6" r="6" fill="currentColor" /><path d="M3.4 6.2 5.2 8l3.4-3.6" stroke="#fff" strokeWidth="1.5" fill="none" /></svg>passes QC</span>
    : <span className="status" style={{ color: "var(--serious)" }}><svg width="12" height="12" viewBox="0 0 12 12"><path d="M6 .8 11.4 10.6H.6z" fill="currentColor" /><path d="M6 4.3v3M6 8.8v.2" stroke="#fff" strokeWidth="1.4" /></svg>flagged{reasons[selected] ? ": " + reasons[selected] : ""}</span>;
  const reads = cell && countLayer ? cell.layers[countLayer] : null;

  return (
    <div className="card">
      {fig && <FigureDialog kind="profile" onClose={() => setFig(false)} />}
      <div className="hd">
        <h2>{name ? <span className="mono">{name}</span> : "Cell profile"}</h2>
        <span className="sub">{status}{ct && <span className="ok"> · {String(ct)}</span>}{view === ALLELE && d.hasAS && <span className="ok"> · allele-specific</span>}</span>
        {xZoomed && <button className="zoomchip" onClick={() => set({ zoom: { x: [0, genomeMax], y: zoom!.y } })} title="Show the whole genome">{zoomText()}  ✕</button>}
        <button className="iconbtn" onClick={() => setFig(true)} disabled={selected == null} title="Download a publication figure"><DlIcon />Figure</button>
      </div>
      {pills.length > 0 && (
        <div className="pills">
          {pills.map(([k, label, v, title]) => (
            <div className={"pill" + (failsFilter(k) ? " bad" : "")} key={k} title={title || metricOf(k)?.help}><span className="k">{label}</span><span className="v">{v}</span></div>
          ))}
        </div>
      )}
      <div className="bd">
        {selected == null ? <div className="empty">Click a row of the heatmap, or pick a cell on the left.</div>
          : !cell ? <div className="empty">Loading…</div>
          : <ZoomBrush>
            <ProfilePlot ref={plot} cell={cell} />
            {d.exprLayer && d.totalLayer && cell.layers[d.exprLayer] && (
              <ExpressionPlot expr={cell.layers[d.exprLayer]} call={cell.layers[d.totalLayer] ?? null}
                callName={d.totalLayer === "total_cn" ? "copy number call" : d.totalLayer} />
            )}
            {d.exprLayer && !d.totalLayer && meta.modality === "scATAC" && cell.layers[d.exprLayer] && (
              <ExpressionPlot expr={cell.layers[d.exprLayer]} call={ternaryOf(cell.layers, d.sem)} callName="call" />
            )}
            {view === ALLELE && d.bafLayer && cell.layers[d.bafLayer] && (
              <BafPlot baf={cell.layers[d.bafLayer]} n={d.bafNLayer ? cell.layers[d.bafNLayer] ?? null : null}
                major={cell.layers.major_cn ?? null} minor={cell.layers.minor_cn ?? null} />
            )}
            {reads && <ReadDepthPlot reads={reads} />}
            <div className="brushhint">drag across the plots to zoom · double-click to see the whole genome</div>
          </ZoomBrush>}
      </div>
    </div>
  );
}

/* ── Lorenz curve of the reads over bins ── */
export function LorenzCard() {
  const { meta, vals, countLayer } = useDerived();
  const selected = useStore((s) => s.selected);
  const theme = useStore((s) => s.theme);
  const cell = useCellDetail(selected);
  const reads = cell && countLayer ? cell.layers[countLayer] : null;
  const gini = selected != null ? vals.gini?.[selected] : undefined;
  const box = useRef<HTMLDivElement>(null);
  const cv = useRef<HTMLCanvasElement>(null);
  const { w: width, h: height } = useSize(box, 260);
  const points = useMemo(() => {
    if (!reads) return null;
    const v = reads.filter((x): x is number => x != null && x >= 0).sort((a, b) => a - b);
    const total = v.reduce((a, b) => a + b, 0);
    if (!v.length || !total) return null;
    const pts: [number, number][] = [[0, 0]];
    let cum = 0;
    v.forEach((x, i) => { cum += x; pts.push([(i + 1) / v.length, cum / total]); });
    return pts;
  }, [reads]);
  useEffect(() => {
    const c = cv.current;
    if (!c || !points) return;
    const g = fit(c, width, height);
    g.clearRect(0, 0, width, height);
    const b = { l: 52, t: 14, w: width - 66, h: height - 58 };
    const ticks = [0, 0.25, 0.5, 0.75, 1], pct = (v: number) => Math.round(v * 100) + "%";
    const { tx, ty } = axes(g, b, ticks, ticks, { x0: 0, x1: 1, y0: 0, y1: 1, xlab: "share of bins, fewest reads first", ylab: "share of reads", xfmt: pct, yfmt: pct });
    g.strokeStyle = css("--ink3"); g.setLineDash([4, 4]); g.beginPath(); g.moveTo(tx(0), ty(0)); g.lineTo(tx(1), ty(1)); g.stroke(); g.setLineDash([]);
    g.fillStyle = css("--accent-soft"); g.beginPath(); g.moveTo(tx(0), ty(0));
    for (const [x, y] of points) g.lineTo(tx(x), ty(y));
    g.lineTo(tx(1), ty(1)); g.closePath(); g.fill();
    g.strokeStyle = css("--accent"); g.lineWidth = 2; g.lineJoin = "round"; g.beginPath();
    points.forEach(([x, y], i) => (i ? g.lineTo(tx(x), ty(y)) : g.moveTo(tx(x), ty(y))));
    g.stroke();
    if (gini != null && !Number.isNaN(gini)) { g.fillStyle = css("--ink"); g.font = `600 12px ${FONT}`; g.textAlign = "left"; g.textBaseline = "top"; g.fillText(`Gini ${gini.toFixed(3)}`, b.l + 10, b.t + 6); }
  }, [points, width, height, gini, theme]);
  void meta;
  return (
    <div className="card">
      <div className="hd"><h2>Coverage evenness</h2><span className="sub">Lorenz curve of the reads over bins</span></div>
      <div className="bd">
        <div ref={box} className="fill" style={{ minHeight: points ? 260 : 0 }}>
          {points ? <canvas ref={cv} /> : <div className="empty">No per-bin read counts in this dataset.</div>}
        </div>
      </div>
    </div>
  );
}

/* ── this cell's QC against all cells ── */
export function CellQCCard() {
  const { metAll, vals, filters } = useDerived();
  const selected = useStore((s) => s.selected);
  const valOf = (k: string) => (selected != null && vals[k] ? vals[k][selected] : NaN);
  return (
    <div className="card">
      <div className="hd"><h2>This cell's QC</h2><span className="sub">against all cells · ▲ beyond median ± 3 MAD on the bad side</span></div>
      <div className="bd tablewrap">
        {metAll.length ? (
          <table>
            <thead><tr><th>metric</th><th>cell</th><th>median</th><th /></tr></thead>
            <tbody>
              {metAll.map((m) => {
                const v = valOf(m.key), out = isOutlier(m, v), fail = fails(filters[m.key], v);
                return (
                  <tr key={m.key}>
                    <td title={m.help}>{m.label}</td>
                    <td style={{ fontWeight: 600, color: fail ? "var(--serious)" : undefined }}>{fmtMetric(m.key, v)}</td>
                    <td className="ok">{fmtMetric(m.key, m.median)}</td>
                    <td style={{ color: "var(--serious)", width: 1, whiteSpace: "nowrap" }}>{fail ? "filter" : out ? "▲" : ""}</td>
                  </tr>
                );
              })}
            </tbody>
          </table>
        ) : <div className="empty">This dataset has no per-cell QC metrics.</div>}
      </div>
    </div>
  );
}
