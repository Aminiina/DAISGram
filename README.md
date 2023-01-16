# DAISGram

DAISGram - C++ Project

Requires Linux OS or Ubuntu VM on Windows
Compiler: g++

## BMP library

In order to test that the (libbmp) library is running correctly,
you need to execute the make testbmp command, then ./test_bmplib
If running correctly, it will generate a checkboard.bmp file.

# Expected results

The 'images' folder contains some images to be used in testing your 
implementation. 
The 'results' folder contains the expected results.

# Example Commands

Help

./main show_help

Gray 

./main images/flower_hires.bmp images/flower_hires.bmp gray out.bmp

Brighten

./main images/fullmoon.bmp images/fullmoon.bmp brighten out.bmp 100


Blend

./main images/blend/blend_a.bmp images/blend/blend_b.bmp blend out.
bmp 0.5 0.5

Concat

./main_tensor tensors/t_5_5_2_progressive.txt tensors/t_5_5_2_random.txt concat out.txt




