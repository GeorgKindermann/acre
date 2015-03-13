# acre
Area-split proportional to pointinfluence 

   Acre was written to calculates growing-space or stand-density of trees
   in a forest. It can also be used to calculate weighted Voroni-Diagrams.
   This means, that the area can be spilt up to the neares point or to the
   point with the greatest influence at a certain place. You can get the
   result as text and visualized (pic format is pbm, pgm ppm).

   The program is written in ANSI-C and is available as the source code so
   you have to compile it before use.

   In Unix you can make this by typing:

   make

   or:

   cc acre.c -lm -o acre

   or:

   gcc acre.c -lm -O2 -o acre

   This produces a program with name acre. The file data.txt is a
   data-sample from a forest. Some help is available by starting the
   program:

   ./acre

   ./acre -h

   ./acre --help

   ./acre -V

   or by reading the source. For copyright-information read also the
   source-code.
