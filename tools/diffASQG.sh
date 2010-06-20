#! /bin/bash
echo $1
diff -y <(zcat $1 | sort) <(zcat $2 | sort) | grep -A2 -B2 ">\|<"
