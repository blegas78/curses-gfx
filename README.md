# curses-gfx
A 3D software renderer using ncurses for terminal graphics

I decided to reinvent the wheel and accru tech debt at the same time by making a software renderer for ascii art.  It's a fun project, almost as challenging and pretty as it useless.

Cmake is used to build installable shared objects.  Recommended install process on something like a Raspberry Pi:

```
mkdir build
cd build
cmake ..
make
sudo make install
```

