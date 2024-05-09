## Test environments

* local OS X install, R 4.2.2
* Windows ( devel, <https://win-builder.r-project.org/> ), R 2023-07-05 r84643 ucrt
* Windows ( devel, windows-x86_64-devel on R-hub )
* Ubuntu Linux ( devel, ubuntu-gcc-devel on R-hub )
* Fedora Linux ( devel, fedora-gcc-devel on R-hub )
* Debian Linux ( devel, debian-gcc-devel on R-hub )


## R CMD check results

0 errors | 0 warnings | 0 note

This result is valid for local OS X install, win-builder devel, and Debian devel

There were 3 notes generated on other platforms (no warnings or errors):

1. 1 warning on both Ubuntu and Fedora devel (R-hub)

```
* checking HTML version of manual ... NOTE
Skipping checking HTML validation: no command 'tidy' found
```

This is a recurring issue on R-hub [R-hub issue #560](https://github.com/r-hub/rhub/issues/548) and can likely be ignored.

2. 2 warnings on the R-hub Windows devel

```
* checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'
```

[R-hub issue #503](https://github.com/r-hub/rhub/issues/503): this seems to be a MiKTeX problem and can likely be ignored.

and

```
* checking for non-standard things in the check directory ... NOTE
Found the following files/directories:
  ''NULL''
```

[R-hub issue #560](https://github.com/r-hub/rhub/issues/560): this seems to be an R-hub issue and so can likely be ignored.