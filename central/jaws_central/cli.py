import click
from jaws_central import config

@click.group()
def jaws():
    pass

@jaws.command()
def start():
    connex = config.connex
    connex.add_api('swagger.yml')

    if __name__ == "__main__":
        connex.run(host='0.0.0.0', port=5000, debug=False)
