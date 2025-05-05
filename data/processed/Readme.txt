brands is the same as data in the raw file.

OP for Rno3273 is based on OP for R of the raw file. 

Compared with OP for R in the raw file, 

OP for Rno3273 has no column for COP. In COP all values are 1. COP is a selection variable. That is, only press releases having COP were included in the dataset.
Due to this selection criteria, we remove the variable COP. 

After removing COP, Documents 9, 21, 33, 70, and 79 become 0 vectors and thus Documents 9, 21, 33, 70, and 79 are removed. CA does not work on a table where 0 vector in the row or column.

In addition, Documents 32 and 73 are removed because they may become 0 vectors if we take (32, Resp.C.I) and (73, PH) as outlying cells.

Finally, we obtain OP for Rno3273