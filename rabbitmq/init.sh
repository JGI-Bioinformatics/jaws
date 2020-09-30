#!/bin/bash

# Create Default RabbitMQ setup
( sleep 10; \
# Create user
rabbitmqctl add_user jaws jawstest; \
# Set user rights
rabbitmqctl set_user_tags jaws administrator; \
# Create vhosts
rabbitmqctl add_vhost jaws_dev;\
rabbitmqctl add_vhost jtm_dev; \
# Set vhost permissions
rabbitmqctl set_permissions -p jtm_dev jaws ".*" ".*" ".*" 
rabbitmqctl set_permissions -p jaws_dev jaws ".*" ".*" ".*" 
) &  
rabbitmq-server $@
