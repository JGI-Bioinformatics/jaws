# To setup dev environment

## Host name in local computer

Globus oauth is configured to work with the dev server at https://jawsd.jgi.doe.gov:3003. So we need to config the server to run with localhost -> jawsd.jgi.doe.gov.

#### Add jawsd.jgi.doe.gov to localhost list

Added line <code> 127.0.0.1 jawsd.jgi.doe.gov </code> to your hosts file (/etc/hosts for Mac and linux). You need sudo for editing this file.

## Create self-signed SSL certificate
Reference [here](https://ksearch.wordpress.com/2017/08/22/generate-and-import-a-self-signed-ssl-certificate-on-mac-osx-sierra/).

Make sure you have <strong>openssl</strong> command installed:
<code>which openssl</code>. If not installed, run <code> brew install openssl</code> (Mac).


Create <strong>server.key</strong>, <strong>server.crt</strong> and <strong>v3.ext</strong> files in <strong>$HOME/.ssh/ssl</strong> directory. 
<br />
Create the .ssh/ssl directory if it is not exits, and move into the directory.

#### server.key

<ul>
  <li> <code>openssl genrsa -des3 -passout pass:<strong>reactjaws</strong> -out server.pass.key 2048 </code>
      <br />
      This creates a file named <code>server.pass.key</code>.  You can use a different value in place of the <code>reactjaws</code> in [pass:]<strong>reactjaws</strong> but should be no shorter than 4 chars
  </li>
  <li> <code>openssl rsa -passin pass:<strong>reactjaws</strong> -in server.pass.key -out server.key </code>
      <br />
      This creates the <code>server.key</code> file.
  </li>
  <li> <code> rm server.pass.key </code>
      <br />
      Do not need this file anymore.
  </li>
</ul>

#### server.crt
##### server.csr - needed to create server.crt

The below command will ask you for information that would be included in the certificate. Since this is a self-signed certificate, there is no need to provide the 'challenge password' (to leave it blank, press enter).
<br />
<code>openssl req -new -key server.key -out server.csr</code>
<pre>
You are about to be asked to enter information that will be incorporated
into your certificate request.
What you are about to enter is what is called a Distinguished Name or a DN.
There are quite a few fields but you can leave some blank
For some fields there will be a default value,
If you enter '.', the field will be left blank.

Country Name (2 letter code) [AU]: <strong>US</strong>
State or Province Name (full name) [Some-State]: <strong>CA</strong>
Locality Name (eg, city) []: <strong>Berkeley</strong>
Organization Name (eg, company) [Internet Widgits Pty Ltd]: <strong>LBL</strong>
Organizational Unit Name (eg, section) []: <strong>SE</strong>
Common Name (e.g. server FQDN or YOUR name) []: <strong>dev.jgi.doe.gov</strong>
Email Address []: <provide the email address to be included in the certificate signing request>
 
Please enter the following 'extra' attributes
to be sent with your certificate request
A challenge password []:
An optional company name []:
</pre>

##### v3.ext - needed to create server.crt
Create <code>v3.ext</code> file with the following content:
<pre>
authorityKeyIdentifier=keyid,issuer
basicConstraints=CA:FALSE
keyUsage = digitalSignature, nonRepudiation, keyEncipherment, dataEncipherment
subjectAltName = @alt_names
 
[alt_names]
DNS.1 = dev.jgi.doe.gov
</pre>

This file is meant to address the page loading error (connection is not private) for Chrome. But it seems does not work. After you setup the SSL properly as the above and still see the loading error, just type <code>thisisunsafe</code> and that will resolve ths issue.

##### create server.crt
<code>openssl x509 -req -sha256 -extfile v3.ext -days 365 -in server.csr -signkey server.key -out server.crt<code>
<pre>server.key -out server.crt
Signature ok
subject=/C=<country>/ST=<state>/L=<locality>/O=<organization-name>/OU=<organization-unit-name>/CN=<common-name-probably-server-fqdn>/emailAddress=<email-address-provided-while-generating-csr>
Getting Private key
</pre>


## Setup your dev web server to run on https://jawsd.jgi.doe.gov:3003

After you have your SSL certificate created properly, make sure you have the following lines in the repo's <code>.env</code> file:

<pre>
HTTPS=true
HOST=jawsd.jgi.doe.gov
PORT=3003
SSL_CRT_FILE=$HOME/.ssh/ssl/server.crt 
SSL_KEY_FILE=$HOME/.ssh/ssl/server.key 

SKIP_PREFLIGHT_CHECK=true
</pre>

<code>npm start</code>
<br />
Once the dev web server is started and running, your browser should pop to foreground and show your page at https://jawsd.jgi.doe.gov:3003.
<br>
If you see <strong>Your connection is not private</strong> page, type <code>thisisunsafe</code> from your keybaord while the page is in focus should return you to the normal app page.
