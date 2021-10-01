# Orbits

This is a simple test program that solves the two body problem ODE system

It was created to test how a library from the [MESA](https://docs.mesastar.org/)
stellar evolution code can be used in another different program. It uses a
Runge-Kutta method that can be found in the `num` library of the MESA src code

To run this code, first cd into `src/mechanics` and compile the code with `./mk`.
Once the code is compiled, run './orbit' and that's it. It needs MESA installed,
which means that you need to set the $MESA_DIR env variable
