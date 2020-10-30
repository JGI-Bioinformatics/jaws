class Channel:
    def __init__(self): pass
    def exchange_declare(self): pass
    def queue_declare(self, queue): pass
    def queue_bind(self): pass
    def basic_qos(self): pass
    def basic_consume(self, queue, auto_ack, on_message_callback): pass
    def start_consuming(self): pass
    def basic_publish(self, exchange, routing_key, body): pass
    def basic_ack(self): pass


class Connection:
    def __init__(self): pass
    def channel(self): return Channel()


class spec:
    def __init__(self): pass

    class Basic:
        def __init__(self): pass
        def Deliver(): pass

    class BasicProperties:
        def __init__(self): pass
