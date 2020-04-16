#!/bin/bash
ls | grep mol2 > /tmp/dock_summarize.$$
while IFS= read -r line
do
  grep --with-filename "Grid_Score:" $line
done < /tmp/dock_summarize.$$
rm /tmp/dock_summarize.$$
