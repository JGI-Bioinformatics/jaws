import axios from "axios"; // this requires our server on this domain

import {cookieGet, cookieSet, cookieRemove} from "./cookie";

export const SSO_BASE = "https://signon.jgi.doe.gov";
export const SSO_URL = `${SSO_BASE}/signon`;
export const DOMAIN = ".jgi.doe.gov";
export const SSO_COOKIE_KEY = "jgi_session";

// axios
export const httpclient = axios.create(
  {
    baseURL: SSO_BASE,
    // withCredentials: true,
    headers: {
      'Content-Type': 'application/json',
      'Accept': 'application/json',
    }
  });

export const clearSSOCookie = () => {
  cookieRemove(SSO_COOKIE_KEY);
};

// so returns to this site after sso logout
export const setSSOReturnURL = () => {
  cookieSet("jgi_return", window.location.href, {domain: DOMAIN, path:"/"});
};

// for SSO signon there will be a cookie of jgi_session=/api/sessions/HASH; return the HASH
export const getSSOHash = () => {
  const sso_session = cookieGet(SSO_COOKIE_KEY);
  if (sso_session){
    const tokens = sso_session.split("/api/sessions/");
    if (tokens.length === 2){
      return tokens[1];
    }
  }
  return null;
}

export const handleLogin = () => {
  setSSOReturnURL();
  window.location.href = SSO_URL;
};