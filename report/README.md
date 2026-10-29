# Complete LaTeX report

`main.tex` is the complete report source. Only **Organization of Work** is left
as a placeholder; enter the group's actual roles and contribution percentages.

The numerical-method discussion follows Lectures 10--13: RK4 stage updates,
Taylor expansions, timestep-normalized truncation error, and absolute stability.
It does not use Butcher conditions, multilinear derivative expansions, or a
Lipschitz error bound. The scalar Taylor calculation illustrates order matching;
the nonlinear model's observed order is demonstrated by the convergence tables.
Aerodynamics, convective heating, radiation, and all numerical results are retained.

## Compile

Upload `main.tex`, `new-aiaa.cls`, and the entire `figures/` folder to the same
Overleaf project and select `main.tex` as the main document. The bibliography
is embedded in `main.tex`; no separate `.bib` file is needed.

With a local TeX installation, run twice from this folder:

```powershell
pdflatex main.tex
pdflatex main.tex
```

Alternatively, `tectonic main.tex` handles its own reruns.
The local build was checked with Tectonic 0.17.0.

## Regenerate the report figures

From the repository root, run:

```powershell
python -m analysis.run_all
python -m analysis.report_figures
```

The first command regenerates numerical outputs. The second reads the CSV
outputs and creates figures with text sized for the report. It does not change
the numerical model or repeat integration. Recompile after regenerating figures.

## Template provenance

The report retains the supplied `new-aiaa` conference format. The bundled class
identifies itself as `new-aiaa.cls`, version 1.2, developed by Overleaf (2018-01-10).
It was obtained from this mirror:
https://raw.githubusercontent.com/ethanleboeuf/AA290/master/new-aiaa.cls
The class file is unmodified.

The code repository URL in the report is the project's existing remote. This
session has not pushed or published changes. `AE370-submission.zip` includes the
local code, numerical outputs, and completed report snapshot for sharing.
