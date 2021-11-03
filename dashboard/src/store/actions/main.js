import * as actionTypes from "./actions"

import { cookieGet, cookieSet, CK_AVATAR_ID } from "../../utils/cookie";

// action creators
export const updateUserAvatarId = (aid) => {
  return {
    type: actionTypes.UPDATE_USER_AVATAR_ID,
    value: aid
  }
}

export const userAvatar = (avatarId) => {
  return dispatch => {
    let aid = cookieGet(CK_AVATAR_ID);
    if (!aid) {
      cookieSet(CK_AVATAR_ID, avatarId);
    } else {
      aid = parseInt(aid);
      dispatch(updateUserAvatarId(aid));
    }
  }
}