
import Cookies from "universal-cookie";
/**
options (object): Support all the cookie options from RFC 6265
   - path (string): cookie path, use / as the path if you want your cookie to be accessible on all pages
   - expires (Date): absolute expiration date for the cookie
   - maxAge (number): relative max age of the cookie from when the client receives it in seconds
   - domain (string): domain for the cookie (sub.domain.com or .allsubdomains.com)
   - secure (boolean): Is only accessible through HTTPS?
   - httpOnly (boolean): Is only the server can access the cookie?
   - sameSite (boolean|none|lax|strict): Strict or Lax enforcement
 */

export const CK_AVATAR_ID = "avatar_id";            // avatar image index

const api = new Cookies();

const defaultOpts = { path: "/" };
export const cookieGet = (key) => api.get(key);
export const cookieSet = (key, value, opts=defaultOpts) => api.set(key, value, opts);
export const cookieRemove = (key, opts=defaultOpts) => api.remove(key, opts);

