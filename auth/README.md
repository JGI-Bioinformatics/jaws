# JAWS Auth

This is an OAuth2 server which provides authentication and user management for all JAWS services.  Records are stored in a relational database.

## Globus

All JAWS users must have Globus accounts as Globus provides the file transfer capabilities.  Not all Globus users are JAWS users; the Globus groups service is used to verify membership in the "jaws_users" group.

## Scopes

Scopes are determine by membership in Globus Groups.  All jaws users belong to the "jaws_users" group and administrators belong to the "jaws_admins" group.
