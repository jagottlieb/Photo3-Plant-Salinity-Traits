"""Write tables into the literature workbook through Excel itself.

openpyxl round-trips drop chart styling and cached formula values in this workbook, so sheets are
added with Excel over COM (Windows + Excel + pywin32 required), which leaves everything else intact.
Excel works on a temporary copy outside OneDrive (AutoSave would otherwise save partial edits),
which replaces the workbook only after the sheet checks pass.
"""
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

XLSX = Path(__file__).resolve().parents[1] / "Photosynthesis Stomatal Conductance Reduction and Dry Weight Data.xlsx"


def _cell(v):
    if isinstance(v, (np.bool_, bool)):
        return bool(v)
    if isinstance(v, (np.integer, np.floating, float, int)):
        return None if pd.isna(v) else float(v)
    return None if v is None else str(v)


def write_sheets(sheets: list[tuple[str, list[str], pd.DataFrame]], xlsx: Path = XLSX) -> None:
    """Replace or append each (name, note lines, table) sheet at the end; notes go above the header."""
    import win32com.client
    tmp = Path(tempfile.mkdtemp()) / xlsx.name
    shutil.copy2(xlsx, tmp)
    excel = win32com.client.gencache.EnsureDispatch(win32com.client.DispatchEx("Excel.Application"))
    excel.Visible = False
    excel.DisplayAlerts = False
    try:
        wb = excel.Workbooks.Open(str(tmp))
        before = [s.Name for s in wb.Sheets]
        new = [name for name, _, _ in sheets]
        for name, notes, df in sheets:
            if name in [s.Name for s in wb.Sheets]:
                wb.Sheets(name).Delete()
            ws = wb.Worksheets.Add(Before=None, After=wb.Sheets(wb.Sheets.Count))
            ws.Name = name
            assert ws.Index == wb.Sheets.Count, f"'{name}' was not added at the end"
            ncol = df.shape[1]
            rows = [[n] + [None] * (ncol - 1) for n in notes]
            rows.append(list(df.columns))
            rows += [[_cell(v) for v in r] for r in df.itertuples(index=False)]
            ws.Range(ws.Cells(1, 1), ws.Cells(len(rows), ncol)).Value = [tuple(r) for r in rows]
            print(f"Wrote {len(df)} rows to '{name}'")
        after = [s.Name for s in wb.Sheets]
        kept = [n for n in before if n not in new]
        assert [n for n in after if n not in new] == kept, f"existing sheets changed: {before} -> {after}"
        wb.Save()
        wb.Close(False)
    finally:
        excel.Quit()
    shutil.copy2(tmp, xlsx)
    shutil.rmtree(tmp.parent, ignore_errors=True)
