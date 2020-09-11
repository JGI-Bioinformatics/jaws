def kill(job_id):
    job_id=params['job_id']
    pid=lookup_pid(job_id)

    list_descendants ()
    {
      local children=$(ps -o pid= --ppid "$1")

      for pid in $children
      do
        list_descendants "$pid"
      done

      echo "$children"
    }

    kill $(list_descendants $$)
