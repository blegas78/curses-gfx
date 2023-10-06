# curses-gfx
A 3D software renderer using ncurses for terminal graphics

This project was originally inspired by a 4D hypercube rotation written in C from [here](https://gist.github.com/Mashpoe/3d949824be514c43b58706eb29c33c43).

I decided to reinvent the wheel and accru tech debt at the same time by making a software renderer for ascii art.  It's a fun project, almost as challenging and pretty as it useless.

Cmake is used to build installable shared objects.  Recommended install process on something like a Raspberry Pi:

```
mkdir build
cd build
cmake ..
make
sudo make install
```


## Examples:

#### cube2
Rasterization using different shaders and three moving light sources

![cube2](https://github.com/blegas78/curses-gfx/blob/main/docs/images/cube2.png?raw=true)

#### stopwatch

![Stopwatch](https://github.com/blegas78/curses-gfx/blob/main/docs/images/stopwatch.png?raw=true)


#### fileLoader
Uses libassimp to load 3D models
![fileLoader](https://github.com/blegas78/curses-gfx/blob/main/docs/images/fileLoader.png?raw=true)

