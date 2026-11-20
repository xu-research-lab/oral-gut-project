# ============================================================
# Shared group configuration for the Oral-gut-cancer mouse
# metagenomic validation project.
# Group key (user-confirmed 2026-08-03, B/D swap):
#   A = Control, D = PD-1, G = PD-1 + E. coli, C = PD-1 + VP,
#   B = PD-1 + Butyrate, E = PD-1 + VP + Butyrate
# Plot order: Control, PD-1, PD-1 + E. coli, PD-1 + VP,
#             PD-1 + Butyrate, PD-1 + VP + Butyrate
# ============================================================
group_order <- c("Control", "PD-1", "PD-1 + E. coli", "PD-1 + VP",
                 "PD-1 + Butyrate", "PD-1 + VP + Butyrate")

group_colors <- c(
  "Control"              = "#E64B35",
  "PD-1"                 = "#4DBBD5",
  "PD-1 + E. coli"       = "#8491B4",
  "PD-1 + VP"            = "#00A087",
  "PD-1 + Butyrate"      = "#3C5488",
  "PD-1 + VP + Butyrate" = "#F39B7F"
)
