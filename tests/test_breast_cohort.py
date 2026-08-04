from __future__ import annotations

from pathlib import Path

import pandas as pd

from reach_gap.breast_cohort import audit_breast_cell_groups


def test_breast_cell_group_audit(tmp_path: Path) -> None:
    first = tmp_path / "mid.csv"
    second = tmp_path / "bottom.csv"
    pd.DataFrame(
        {"cell_id": ["a", "b", "c"], "group": ["0_Tumor_Cells", "1_T_Cells", "2_Tumor_Cells"]}
    ).to_csv(first, index=False)
    pd.DataFrame({"cell_id": ["d", "e"], "group": ["0_Stromal_Cells", "1_Tumor_Cells"]}).to_csv(
        second, index=False
    )
    result = audit_breast_cell_groups({"mid": first, "bottom": second}, tmp_path / "out")
    assert result["cells"] == 5
    assert result["spatial_expression"]["status"] == "NOT_COMPUTED"
    composition = pd.read_csv(tmp_path / "out" / "breast_cell_group_composition.csv")
    assert "Tumor_Cells" in set(composition["broad_group"])
