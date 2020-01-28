# worker

JAWS Worker node runs alongside Cromwell and JTM Manager on a computing resource (cluster).  The JAWS Central server communicates with the workers via RabbitMQ.  The Worker manages the transfer of files (via Globus) and the running of workflows (via Cromwell).  This package provides no user authentication, which is handled upstream by JAWS Central.
