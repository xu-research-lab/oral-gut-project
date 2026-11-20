# ============================================================
# Shared group definitions (user-confirmed 2026-08-02, rev 2)
# Labels: PD-1 (with hyphen). Order:
#   Control -> PD-1 -> PD-1 + E. coli -> PD-1 + VP
#   -> PD-1 + VP + Butyrate -> PD-1 + Butyrate (last)
# ============================================================
group_map <- c(
  VP1 = "PD-1 + VP", VP2 = "PD-1 + VP", VP3 = "PD-1 + VP", VP4 = "PD-1 + VP", VP5 = "PD-1 + VP",
  E1  = "PD-1 + E. coli", E2 = "PD-1 + E. coli", E3 = "PD-1 + E. coli", E4 = "PD-1 + E. coli", E5 = "PD-1 + E. coli",
  K1  = "Control", K2 = "Control", K3 = "Control", K4 = "Control", K5 = "Control",
  P1  = "PD-1", P2 = "PD-1", P3 = "PD-1", P4 = "PD-1", P5 = "PD-1",
  H1  = "PD-1 + VP + Butyrate", H2 = "PD-1 + VP + Butyrate", H3 = "PD-1 + VP + Butyrate",
  H4  = "PD-1 + VP + Butyrate", H5 = "PD-1 + VP + Butyrate",
  J1  = "PD-1 + Butyrate", J2 = "PD-1 + Butyrate", J3 = "PD-1 + Butyrate",
  J4  = "PD-1 + Butyrate", J5 = "PD-1 + Butyrate"
)
group_order <- c("Control", "PD-1", "PD-1 + E. coli", "PD-1 + VP", "PD-1 + VP + Butyrate", "PD-1 + Butyrate")
group_colors <- c(
  "Control" = "#4DAF4A",
  "PD-1" = "#984EA3",
  "PD-1 + E. coli" = "#377EB8",
  "PD-1 + VP" = "#E41A1C",
  "PD-1 + VP + Butyrate" = "#FF7F00",
  "PD-1 + Butyrate" = "#A65628"
)
