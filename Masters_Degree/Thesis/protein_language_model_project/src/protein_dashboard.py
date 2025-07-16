import pandas as pd
import sqlite3
from shiny import App, render, ui
import sys

# -------------------------------------------------
# Smart platform-based database path
# -------------------------------------------------
if sys.platform.startswith("linux"):
    # WSL or native Linux
    db_path = "/mnt/d/protein_hallucination_db/protein_hallucination.db"
else:
    # Windows native (Spyder / Anaconda)
    db_path = r"D:\protein_hallucination_db\protein_hallucination.db"

# -------------------------------------------------
# Load data
# -------------------------------------------------
def load_data():
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query("SELECT * FROM results", conn)
    conn.close()
    return df

df = load_data()

# -------------------------------------------------
# Prepare family list for dropdown
# -------------------------------------------------
family_list = df['family'].dropna().unique().tolist()
pdb_list = df['pdb_id'].dropna().unique().tolist()
# -------------------------------------------------
# Shiny UI
# -------------------------------------------------
app_ui = ui.page_fluid(
    ui.tags.style("""
        table.dataframe {
            font-family: "Courier New", Courier, monospace;
        }
    """),
    ui.h2("Protein Language Model Dashboard"),
    ui.input_select(
        "family_filter",
        "Select Family",
        choices=family_list,
        selected=family_list[0] if family_list else None
    ),
    ui.input_select("pdb_id_filter",
                    "select_pdb_id",
                    choices=pdb_list,
                    selected=pdb_list[0] if pdb_list else None),
    ui.row(
        ui.column(4,
            ui.card(
                ui.h4("Sequence"),
                ui.output_text("sequence_text")
            )
        ),
        ui.column(8,
            ui.output_table("filtered_table")    
        )
    )
    
)

# -------------------------------------------------
# Shiny server
# -------------------------------------------------
def server(input, output, session):
    @output()
    @render.table
    def filtered_table():
        selected_family = input.family_filter()
        selected_pdb_id = input.pdb_id_filter()
        filtered = df[(df['family'] == selected_family) & (df['pdb_id'] == selected_pdb_id)]
        if not filtered.empty:
            return filtered[['pdb_id', 'wt_or_var','sequence' ]]
        else:
            return pd.DataFrame(columns=['pdb_id', 'wt_or_var','sequence'])
    @output()
    @render.text
    def sequence_text():
        selected_family = input.family_filter()
        filtered = df[df['family'] == selected_family]
        if not filtered.empty:
            return str(filtered.iloc[0]['sequence'])
        else:
            return "No sequence found."

# -------------------------------------------------
# Create app
# -------------------------------------------------
app = App(app_ui, server)
