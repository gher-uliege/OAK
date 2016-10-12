
while true; do
  change=$(inotifywait -e close_write,moved_to,create *.F90 *.f90)

  if ./test_toymodel; then
     echo OK
     notify-send -t 1000 -u low 'OK' 'Compilation OK'
  else
     echo FAIL
     notify-send -t 2000 -u low 'Fail' 'Compilation fails'
  fi
done
