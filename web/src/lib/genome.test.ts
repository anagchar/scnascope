import { describe, expect, it } from "vitest";
import { barEnd, binColumns, typicalSpan } from "./genome";

// chr1: bins 0-10, 10-20, (gap 20-40), 40-50; chr2 after a gap: 60-70
const B = { chr: ["1", "1", "1", "2"], start_cum: [0, 10, 40, 60], end_cum: [10, 20, 50, 70] };

describe("binColumns", () => {
  it("fills a gap inside a chromosome with the nearer bin", () => {
    // pixel centres at 5, 15, ..., 75
    expect(Array.from(binColumns(B, 0, 80, 8))).toEqual([0, 1, 1, 2, 2, -1, 3, -1]);
  });
});

describe("barEnd", () => {
  it("draws a bin far wider than typical at a typical width, others whole", () => {
    const S = [0, 10, 20, 30], E = [10, 20, 30, 330]; // the last bin spans 300
    const t = typicalSpan(S, E);
    expect(t).toBe(10);
    expect(barEnd(S, E, 1, t)).toBe(20);
    expect(barEnd(S, E, 3, t)).toBe(40);
  });
});
