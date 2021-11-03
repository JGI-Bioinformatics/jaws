import * as actionTypes from "./actions";

import { getUser } from "../../utils/HttpClientProvider";
import { getSSOHash } from "../../utils/sso";

export const authStart = () => {
  return {
      type: actionTypes.AUTH_START
  }
}

export const authSuccess = (user, kcUser) => {
  // console.log(user, '[auth.js] authSuccess')
  const access_token = user.jaws_token;
  return {
      type: actionTypes.AUTH_SUCCESS,
      user: {
        first_name : kcUser ? kcUser.given_name : user.sso_json.user.first_name,
        last_name : kcUser ? kcUser.family_name : user.sso_json.user.last_name,
        email_address : kcUser ? kcUser.email : user.sso_json.user.email_address,
        // email : kcUser ? kcUser.email : user.sso_json.user.email_address,
      },
      access_token
  }
}

export const authFail = (error) => {
  return {
      type: actionTypes.AUTH_FAIL,
      error: error
  }
}

export const logout = () => {
  return {
      type: actionTypes.AUTH_LOGOUT
  }
}

/*
  async:
  sso_hash: the token value from SSO cookie value : /api/sessions/TOKEN
  gkUser: keycloak user; Layout detects and sets value when calling auth()
*/
export const auth = (kc=null) => {
  let authVal = getSSOHash();
  let kcUser = null;

  if(kc && kc.authenticated){
    // console.log(kc.tokenParsed, "[auth] kc");
    authVal = kc.tokenParsed.email;
    kcUser = kc.tokenParsed;
  }
  return dispatch => {
    
    // console.log(authVal, "[auth]");
    if (authVal){
      dispatch(authStart())
      getUser(authVal)
      .then((user) => {
        dispatch(authSuccess(user, kcUser))
      })
      .catch((error) => {
        dispatch(authFail(error))
      })
    }
  }
}