# JAWS Auth

This is an OAuth2 server which provides authentication and user management for all JAWS services.  Records are stored in a relational database.

## Globus

All JAWS users must have Globus accounts as Globus provides the file transfer capabilities.  Obviously, not all Globus users are JAWS users.  This service automatically authenticates against the Globus OAuth2 service, which, in turn, may authenticate against any one of several supported OAuth2 services (e.g. Google, OrcID, etc.).

## Scopes

In addition to the user information provided by Globus, JAWS Auth also stores the JAWS user's access level.

- user : a registered user is able to submit analysis runs and manage workflows in the catalog which they own
- moderator : is effectively an owner of all workflows in the catalog; this allows staff to fix a broken WDL, for example.
- admin : has access to all user's analyses and stats; for example, an admin may manage cancel any user's analysis or change the owner of an analysis
