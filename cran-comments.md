# Resubmission: fixed the DESCRIPTION citation and reduced the example run time

## Test environments
* local Windows 11, R 4.5.2
*
  win-builder (devel and release)
* macOS
  (mac-builder)
* R-hub: ubuntu, ...

## R CMD check
  results
0 errors  | 0 warnings  | 0 notes 
* New
  submission.

## Notes
* Examples that need
  internet or run a long time are wrapped in
  \donttest{}.
* Parallel code uses at most 2 cores
  when running checks.
