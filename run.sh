#!/bin/bash

echo "Run handin template"

echo "Creating the plotting directory if it does not exist"
if [ ! -d "plots" ]; then
  echo "Directory does not exist, creating it!"
  mkdir plots
fi

# Run poisson.py script
echo "Run the Satellite script ..."
python3 satellite.py

# Run vandermonde.py script
echo "Run the Heating script ..."
python3 heating.py

# Run your other Python scripts as needed
# ...


echo "Generating the PDF"

pdflatex hand-in-2.tex
bibtex hand-in-2.aux
pdflatex hand-in-2.tex
pdflatex hand-in-2.tex


echo "Script execution completed."
