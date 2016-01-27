TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt


QMAKE_CFLAGS += -std=c99 -funroll-all-loops -ftree-vectorizer-verbose=10 --param max-unroll-times=32
SOURCES += main.c

