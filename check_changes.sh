
while true; do
  change=$(inotifywait -e close_write,moved_to,create .)

  if ./toymodel_check; then
     echo OK
     notify-send -t 1000 -u low 'OK' 'Compilation OK'
  else
     echo FAIL
     notify-send -t 2000 -u critical 'Fail' 'Compilation fails'
  fi
done
