#!/usr/bin/python

# Print statement that behaves the same for python 2 and 3.
# E,g, print_(1.0, 2, "5") will always print the string "1.0 2 5".
def print_(*args):
    print(" ".join(str(a) for a in args))
