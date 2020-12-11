# AMQPStorm-RPC

A JSON-RPC2 library based on AMQPStorm.

## Includes:

*  rpc_server.py : used to process JSON-RPC2 requests
*  rpc_client.py : used to send JSON-RPC2 requests to the rpc_server
*  rpc_index.py  : used to manage a collection of rpc clients.

### RPC Server:

#### Synopsis:


#### Description:


### RPC Client:

#### Synopsis:


#### Description:


### RPC Index:

#### Synopsis:

`  import amqpstorm_rpc.rpc_index
`  import configparser
`
`  config_file = './rpc_clients.ini'
`  config = configparser.ConfigParser()
`  config.read(config_file)
`  rpc_index = amqpstorm_rpc.RPC_Index(config)
`  client = rpc_index.get_client('EAGLE')
`  try:
`    response = client.request("hello", { "name": "Larry" })
`  except:
`    raise
`  if "error" in response:
`    raise SystemExit(f'Error {(response["error"]["code"]}: {response["error"]["message"]}')
`  print(response["result"])

#### Description:

This class is only required if you have several clients to keep track of.  Clients may use different AMQP services or have unique credentials.



## Credits:

A JSON-RPC2 library based upon AMQPStorm by E. Andersson.

AMQPStorm was written by E. Andersson and this module is largely his RPC example with minor extensions by the JAWS Team.

See:
  - https://www.amqpstorm.io/
  - https://github.com/eandersson/amqpstorm/tree/master/examples
