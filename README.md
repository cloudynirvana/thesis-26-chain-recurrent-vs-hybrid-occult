# Hybrid occult mode switches versus chain-recurrent components on a joint hybrid–smooth field

**Thesis #26.** Computational research, set out in Nile University B.Sc. chapter order for handoff.

**Depends on:** Thesis #4 (hybrid occult modes) and Thesis #17 (complete Lyapunov / chain-recurrent partition).

**Author:** Kelechi Emeka Ogbonna  
**Email:** kelechiogbonna300@gmail.com  
**GitHub:** https://github.com/cloudynirvana  
**Date:** 21 September 2026

On a joint hybrid–smooth toy, do hybrid occult mode-switch observables align with chain-recurrent components from a complete Lyapunov construction, or can switches fall inside a single recurrent region?

They do not align. The smooth proliferative field has an unstable focus and a periodic orbit of period 19.876. Those are two chain-recurrent pieces. An immune guard and an angiogenic guard each meet the orbit twice, at points at least 0.122 from the focus, so a switch can sit in the orbit alone. The same angiogenic line is also crossed in the gap between the focus and the orbit. A cycling guard at T = b misses the orbit (the orbit reaches only T = 0.418) and a transient from (0.55, 0.80) crosses it outside the collocation defect. On 22×22, 34×34 and 42×42 grids the four orbit hits fall in the interior of one defect component under a predeclared radius rule.

The defect is a collocation defect. It is not a certified Conley set. After one rebuild the orbital derivative is still positive on 47.6 percent of the orbit, and the peak-to-peak change of V on that orbit is 0.0243 against a transient drop of 2.109. No number in Chapter Four is copied from the T04 or T17 deposits. Both are recomputed here only in the sense that the mode names and the collocation formulae are theirs; the arithmetic is this script.

This is research only. It is not a medical device, not clinical decision support, not a dose, and not a cure. No document DOI is registered.

See [DISCLAIMER.md](DISCLAIMER.md). The manuscript is [THESIS.md](THESIS.md).

## Files

| Path | Role |
| --- | --- |
| `THESIS.md` | Manuscript (Chapters 1 to 5, Vancouver citations) |
| `THESIS.pdf` | PDF built from the Markdown |
| `build_pdf.py` | Regenerates `THESIS.pdf` |
| `CITATION.cff` | Citation metadata, no document DOI |
| `DISCLAIMER.md` | Research-only boundary |
| `sim/joint.py` | Seed-labelled joint toy (run label 20260921; the arithmetic is deterministic) |
| `sim/results.json` | Numbers cited in Chapter Four |
| `sim/figures/` | Partition, guard classes, orbital derivative, level set, decrease, hybrid traces |

## Reproduce

```bash
python3 -m pip install -r sim/requirements.txt
python3 sim/joint.py
python3 build_pdf.py
```

NumPy, SciPy and Matplotlib are required for the toy. The PDF step also needs the `markdown` and `weasyprint` packages. Regenerating the script rewrites `sim/results.json` and `sim/figures/`.

## Cite

Ogbonna KE. Hybrid occult mode switches versus chain-recurrent components on a joint hybrid–smooth field [Internet]. Thesis #26 computational research thesis. 21 September 2026 [cited YYYY Mon DD]. Available from: https://github.com/cloudynirvana/thesis-26-chain-recurrent-vs-hybrid-occult

Machine-readable fields are in `CITATION.cff`. Add a document DOI there only after one exists.

Hub index, for cataloguing only: [research-theses-hub](https://github.com/cloudynirvana/research-theses-hub).

## Licence

Text and sketch code are MIT, with attribution. Computational research only.
