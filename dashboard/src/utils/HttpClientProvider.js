import axios from "axios";
/*
 * http interface to jaws backend server
 * swagger : http://jaws.lbl.gov:5002/api/v2/ui/#/default
 */

// JAWs backend api server URL
/*
  https://jaws.lbl.gov:8011 -> https://jaws.lbl.gov:5001; dev
  https://jaws.lbl.gov:8012 -> https://jaws.lbl.gov:5002; dev
  https://jaws.lbl.gov:8013 -> https://jaws.lbl.gov:5003; dev
 */
const JAWS_API_SERVER_DEV = "https://jaws.lbl.gov:8011"; // process.env.REACT_APP_BACKEND_API_SERVER || "https://jaws.lbl.gov:8011";
// const JAWS_API_SERVER_STAGING = "http://jaws.lbl.gov:5002";
// const JAWS_API_SERVER_PROD = "http://jaws.lbl.gov:5003";

const API_PREFIX="/api/v2";
export const JAWS_API_STATUS = `${API_PREFIX}/status`;             // the status endpoint
export const JAWS_API_AUTH = `${API_PREFIX}/auth`;
export const JAWS_API_WDL = `${API_PREFIX}/wdl`;
export const JAWS_API_SEARCH = `${API_PREFIX}/search`;
export const JAWS_API_RUN = `${API_PREFIX}/run`;
export const JAWS_API_SITE = `${API_PREFIX}/site`;

// axios
const httpclient = axios.create(
  {
    baseURL: JAWS_API_SERVER_DEV,
    // withCredentials: true,
    headers: {
      'Content-Type': 'application/json',
      'Accept': 'application/json',
    }
  });

export const getUser = (authVal) => {
  const api = `${JAWS_API_AUTH}/${authVal}`;
  const headers = { Authorization: `Bearer ${process.env.REACT_APP_JAWS_API_MASTER}` }
  // console.log(api, process.env.REACT_APP_JAWS_API_MASTER, "[HttpClientP] getUser");
  return httpclient.get(api, { headers: headers })
            .then(resp => {
              /*
              resp.data
              {
                "jaws_token": JAWS_USER_TOKEN,
                "sso_json": {
                  "user": { // FROM SSO service directly
                    "address_1": ADDR,
                    "address_2": "",
                    "city": CITY,
                    "comments": null,
                    "contact_id": #####,
                    "contact_update_required": null,
                    "country": COUNTRY,
                    "created_at": TS,
                    "department": null,
                    "email": EMAIL,
                    "email_address": EMAIL,
                    "failed_login_at": null,
                    "fax_number": null,
                    "first_name": FIRST_NAME,
                    "funding_sources": "",
                    "gender": null,
                    "id": #####,
                    "institution": INS,
                    "institution_title": "",
                    "institution_title_comment": "",
                    "institution_type": INST_TYPE,
                    "internal": true,
                    "kbase_username": "",
                    "last_authenticated_at": TS,
                    "last_name": LAST_NAME,
                    "login": LOGIN_USER_NAME,
                    "middle_name": null,
                    "num_failed_logins": 0,
                    "orcid_id": null,
                    "phone_number": ###-###-####,
                    "postal_code": #####,
                    "prefix": null,
                    "replaced_by": null,
                    "replaces": [],
                    "small_business": null,
                    "state": "CA",
                    "suffix": null,
                    "updated_at": TS
                  }
                }
              }
              */
              return resp.data 
            })
            .catch((error) => {
              return null;
            })
}

export default httpclient;
