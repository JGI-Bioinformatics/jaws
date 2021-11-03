import * as actionTypes from "../actions/actions";
import { updateObject } from "../utility";

const initState = {
  user : null,
  access_token : null, //access token from /v2/oauth2/token API, for jaws endpoint access
  loading: false,
  error: null
}

const authStart = (state) => {
  return updateObject(state, { loading: true })
}

const authSucceed = (state, action) => {
  return updateObject(state, {
      loading: false,
      user: action.user,
      access_token: action.access_token,
      error: false,
      ...action.user
  })
}

const authFail = (state, action) => {
  return updateObject(state, {
      loading: false,
      error: action.error,
  })
}

const authLogout = (state) => {
  return updateObject(state, { 
    user: null, 
    access_token: null, 
    loading: false, 
    error: null 
  })
}

const reducer = (state=initState, action) => {
  switch (action.type){
    case actionTypes.AUTH_START: return authStart(state)
    case actionTypes.AUTH_SUCCESS: return authSucceed(state, action)
    case actionTypes.AUTH_FAIL: return authFail(state, action)
    case actionTypes.AUTH_LOGOUT: return authLogout(state)
    default: return state
  }
}

export default reducer;