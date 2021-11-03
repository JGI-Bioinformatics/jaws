import * as actionTypes from "../actions/actions";
import {cookieSet, CK_AVATAR_ID} from "../../utils/cookie";

const initState = {
  avatarId :0
}

const reducer = (state=initState, action) => {
  switch (action.type){
    case actionTypes.UPDATE_USER_AVATAR_ID:
      cookieSet(CK_AVATAR_ID, action.value);
      return {
        ...state,
        avatarId: action.value
      }
    default:
      return state;
  }
}

export default reducer;