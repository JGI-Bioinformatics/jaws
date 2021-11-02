import React from 'react';
import { ReactKeycloakProvider, useKeycloak } from "@react-keycloak/web";
import ReactDOM from 'react-dom';
import './index.css';
import './App.css';
import App from './App';
import keycloak from "./utils/keycloak";

const AppMain = () => {
  const { initialized } = useKeycloak();

  if (!initialized) {
    return <div style={{ padding: "20px" }}>Loading...</div>;
  }

  return (<App /> );
};

ReactDOM.render(
  <React.StrictMode>
    <ReactKeycloakProvider authClient={keycloak}>
      <AppMain />
    </ReactKeycloakProvider>
  </React.StrictMode>,
  document.getElementById('root')
);

// If you want your app to work offline and load faster, you can change
// unregister() to register() below. Note this comes with some pitfalls.
// Learn more about service workers: https://bit.ly/CRA-PWA
// serviceWorker.unregister();
