# curses-gfx
A 3D software renderer using ncurses for terminal graphics

This project was originally inspired by a 4D hypercube rotation written in C from [here](https://gist.github.com/Mashpoe/3d949824be514c43b58706eb29c33c43).

I decided to reinvent the wheel and accru tech debt at the same time by making a software renderer for ascii art.  It's a fun project, almost as challenging and pretty as it useless.

## Dependencies
```
ncurses is the only required dependency.  
sudo apt install libncurses5-dev    # Linux, like Raspbian
brew install ncurses                # macos with Brew
```

The optional dependencies are really nice to have to expand capability.

```
sudo apt install libassimp-dev libbullet-dev libsdl2-dev    # Linux, like Raspbian
brew install assimp bullet sdl2                             # macos with Brew
```

## Setup
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

#### cubepng
Loads a PNG file to be used as a model texture
![fileLoader](https://github.com/blegas78/curses-gfx/blob/main/docs/images/cubepng.png?raw=true)


#### fileLoader
Uses libassimp to load 3D models
![fileLoader](https://github.com/blegas78/curses-gfx/blob/main/docs/images/fileLoader.png?raw=true)

![fileLoader2](https://github.com/blegas78/curses-gfx/blob/main/docs/images/mario.png?raw=true)

#### bullet
An example integrating the Bullet physics engine.  Watch the following video in the highest resolution setting if possible:

[![bullet](https://img.youtube.com/vi/6goxQYAVHiQ/maxresdefault.jpg)](https://youtu.be/6goxQYAVHiQ)

#### stopwatch

![Stopwatch](https://github.com/blegas78/curses-gfx/blob/main/docs/images/stopwatch.png?raw=true)
