#!/usr/bin/env python3
"""Click at screen coordinates using Quartz (pyobjc).

Usage: /usr/bin/python3 _click_at.py <x> <y>

System Python at /usr/bin/python3 bundles pyobjc (Quartz). Used as a fallback
when System Events cannot invoke Java Swing widgets via AX actions and we need
to click an exact pixel position.
"""
import sys
import time
import Quartz


def click_at(x: float, y: float) -> None:
    pt = Quartz.CGPointMake(x, y)
    # Move
    mv = Quartz.CGEventCreateMouseEvent(None, Quartz.kCGEventMouseMoved, pt, 0)
    Quartz.CGEventPost(Quartz.kCGHIDEventTap, mv)
    time.sleep(0.05)
    # Down
    dn = Quartz.CGEventCreateMouseEvent(
        None, Quartz.kCGEventLeftMouseDown, pt, Quartz.kCGMouseButtonLeft
    )
    Quartz.CGEventPost(Quartz.kCGHIDEventTap, dn)
    time.sleep(0.05)
    # Up
    up = Quartz.CGEventCreateMouseEvent(
        None, Quartz.kCGEventLeftMouseUp, pt, Quartz.kCGMouseButtonLeft
    )
    Quartz.CGEventPost(Quartz.kCGHIDEventTap, up)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: _click_at.py <x> <y>", file=sys.stderr)
        sys.exit(1)
    click_at(float(sys.argv[1]), float(sys.argv[2]))
