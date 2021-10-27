import React from "react";
import { ReactQueryDevtools } from "react-query-devtools";

import { StylesProvider } from "@material-ui/core/styles";
import { BrowserRouter } from "react-router-dom";

import { Provider } from "react-redux";
import { createStore, combineReducers, applyMiddleware, compose } from "redux";
import thunk from "redux-thunk";

import Layout from "./components/Layout/Layout";
import { hcSendSpan, hcBuildTrace, hcRootTraceId } from "./utils/honeycomb";

//- use redux to store our states shared by many components
import mainReducer from "./store/reducers/main";
import authReducer from "./store/reducers/auth";

// import withLD from "./components/hoc/withLD"

const rootReducer = combineReducers(
  {
    main: mainReducer,
    auth: authReducer
  }
)

const composeEnhancers = (process.env.REACT_APP_REDUX_DEV && window.__REDUX_DEVTOOLS_EXTENSION_COMPOSE__) || compose;
const store = createStore(rootReducer, /* preloadedState, */ composeEnhancers(
  applyMiddleware(thunk)
));

function App(props) {
  // console.log(props, "App")
  const payload = {
    name: "Application Startup",
    service_name: "startup",
  };
  const trace = hcBuildTrace(
    hcRootTraceId,
    hcRootTraceId,
    null
  );
  hcSendSpan(payload, trace);

  return (
    <>
      <Provider store={store}>
        <StylesProvider injectFirst> {/* so my css modules have higher specificity over @material */}
          <BrowserRouter>
              <Layout {...props}/>
          </BrowserRouter>
        </StylesProvider>
      </Provider>
      {process.env.REACT_APP_QUERY_DEV && <ReactQueryDevtools initialIsOpen={false}/>}
    </>
  );
}

export default App;
