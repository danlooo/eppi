#!/usr/bin/env python3

#
# centrifuge taxon count to lineage cont
#

import pandas as pd
import sqlite3
import sys

tax_count_file = sys.argv[1]
tax_db_file = sys.argv[2]
out_file = sys.argv[3]

tax_df = pd.read_csv(tax_count_file, header=None)
tax_df.columns = ["hit_count", "tax_id"]
tax_df = tax_df.set_index("tax_id")

cur = sqlite3.connect(tax_db_file).cursor()

for r_id, r in tax_df.iterrows():
	try:
		lineage = cur.execute("select * from lineage where tax_id == '%s' limit 1;" % r_id).fetchone()
		tax_df.loc[r_id, "species"] = lineage[1]
		tax_df.loc[r_id, "genus"] = lineage[2]
		tax_df.loc[r_id, "family"] = lineage[3]
		tax_df.loc[r_id, "order"] = lineage[4]
		tax_df.loc[r_id, "class"] = lineage[5]
		tax_df.loc[r_id, "phylum"] = lineage[6]
		tax_df.loc[r_id, "superkingdom"] = lineage[7]
	except:
		print("taxon %s passed." % r_id)

tax_df.to_csv(out_file)
