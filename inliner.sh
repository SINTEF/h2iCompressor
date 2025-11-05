#!/bin/bash

file="notebook.py"
pkgdir="h2iCompressor"

if [[ -z "$file" || ! -f "$file" ]]; then
  echo "Error: your must specify an exisiting file to process"
  exit 1
fi
if  [[ ! -d $pkgdir ]]; then
  echo "Error: your must specify a directory which has the source files"
  exit 1
fi

linenum=0
declare -a linenums
declare -a inserts

while IFS= read -r line; do
  ((linenum ++))
  if [[ "$line" == *"#inline "* ]]; then
    linenums+=($linenum)
    linearr=($line)
    inserts+=(${linearr[1]})
  fi
done < "$file"

if [[ $1 == "-f" ]]; then
  response="y"
else
  echo "The program will perform the following inlining:"
  echo "Files to insert:" "${inserts[@]}"
  echo "After line numbers:" "${linenums[@]}"

  read -r -p "Do you want to proceed? [y/N] " response
fi

if [[ "$response" =~ ^(yes|y)$ ]]; then
    cp $file $file.bak

    for ins in "${inserts[@]}"; do
      if [[ ! -f $pkgdir/$ins ]]; then
        echo "Error: could not find $pkgdir/$ins, skipping..."
      else
        sed -i -e "/#inline $ins/r $pkgdir/$ins" $file
        sed -i "/#inline $ins/#Have inlined $ins/" $file
      fi
    done

    mv $file output_$file
    mv $file.bak $file
else
    echo "Cancelled the process."
    exit 1 
fi
