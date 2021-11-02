
import React, { useEffect } from "react";
import { useKeycloak } from "@react-keycloak/web";

// import { useHistory, useLocation } from "react-router-dom";

import { Grid } from "@material-ui/core";
import CircularProgress from '@material-ui/core/CircularProgress';
import { Route, withRouter, Switch } from "react-router-dom";

import { useSelector, useDispatch } from 'react-redux';


// import ContentPanel from "../Contents/ContentPanel";
import Box from "@material-ui/core/Box";
import Search from "../Contents/Search/Search";
import SystemInfo from "../Contents/SystemInfo/SystemInfo";
import JobInfo from "../Contents/JobInfo/JobInfo";
import Workflow from "../Contents/Workflow/Workflow";
import User from "../Contents/User/User";
import Advanced from "../Contents/PipelineBuilder/PipelineBuilder";
import MyJobs from "../Contents/JobInfo/MyJobs";

import BlastDB from "../Contents/BlastDB/BlastDB";

import Header from "../Header/Header";
import Footer from "../Footer/Footer";

import Controls from "../Menu/Controls";
import withLD from "../hoc/withLD";

import { setSSOReturnURL } from "../../utils/sso";
import { SSO_URL } from "../../utils/sso";

import { auth as doAuth, logout as onLogout } from "../../store/actions/auth";
import { userAvatar as doAvatarId } from "../../store/actions/main";

import { handleLogin as doLogin } from "../../utils/sso";

import classes from "./Layout.module.css";

const Layout = (props) => {
  
  const { keycloak } = useKeycloak();
  const {keycloak: keycloak_ff, builder: builder_ff} = props.flags;
  // console.log(props.flags, "[Layout]")

  const avatarId = useSelector((state) => state.main.avatarId);
  const loading = useSelector((state) => state.auth.loading);
  const user = useSelector((state) => state.auth.user);
  const dispatch = useDispatch();

  const pageMap = [
    {
      label: "My Jobs", 
      path: "/", 
      component: MyJobs,
      inuse: user ? true : false,
    },
    {
      label: "Site Search", 
      path: user ? "/search" : "/", 
      component: Search,
      inuse: true,
    }, 
    {
      label: "System info", 
      path: "/sysinfo", 
      component: SystemInfo,
      inuse: true,
    }, 
    {
      label: "Job info", 
      path: "/jobinfo", 
      component: JobInfo,
      inuse: false,
    }, 
    {
      label: "Workflow Catalog", 
      path: "/workflow", 
      component: Workflow,
      inuse: true,
    },
    {
      label: "BlastDB Catalog", 
      path: "/blastdb", 
      component: BlastDB,
      inuse: false,
    },
    {
      label: "Pipeline Builder", 
      path:"/builder", 
      component:Advanced,
      inuse: builder_ff ? true : false,
    }
  ]


  const handleLogout = () => {
    if(keycloak_ff){
      keycloak.logout();
    } else {
      setSSOReturnURL();
      window.location.href = `${SSO_URL}/destroy`;
      dispatch(onLogout());
    }
  }

  const handleLogin = () => {
    if(keycloak_ff){
      if (!keycloak.authenticated) {
        keycloak.login();
        dispatch(doAuth(keycloak));
      } 
    } else {
      doLogin();
    }
  }

  const handleUserAccountInfo = () => {
    // honeycomb.sendUiInteractionSpan(`clicked-external-link: edit-jgi-user`);
    if(keycloak_ff){
      keycloak.accountManagement();
    } else {
      window.location.href = "https://signon.jgi.doe.gov/user/contacts/edit";
    }
  }

  // console.log(keycloak, keycloak_ff, user, "[Layout]");

  useEffect(()=>{
    dispatch(doAvatarId(avatarId));
    const kc = keycloak_ff ? keycloak : null;
    dispatch(doAuth(kc));
  }, [avatarId, dispatch, keycloak_ff, keycloak])

  const main = () => {
    // keep only inuse menual items
    const pmap = pageMap.filter( p => p.inuse);
    // console.log(pmap, "[Layout] pmap");

    // console.log(pmap, "[routes] map")
    // need to pass flags in props down to page component ...
    const routes = pmap.map( p => <Route key={p.label} path={p.path} exact 
      render={ () => { return <p.component {...props}/>} } />)

  
    // console.log(routes, "[routes] comp")
    return (
      <>
        <Grid container justify="flex-start">
          <Grid item xl={2}>
            <Controls controls={pmap.map(p=>[p.label, p.path])} {...props} user={user} />
          </Grid>  
          <Grid item xs={10} sm={10} md={9} lg={9} xl={9}>
            <div>
              <Box p={2}>
                <Switch>
                  { routes }
                  {/* <Route>
                    <Redirect to="/" />;
                  </Route> */}
                </Switch>
              </Box>
            </div>
          </Grid>
        </Grid>
      </>
    )
  }

  const content = loading ? <CircularProgress /> :
    (
      <Switch>
        <Route path="/user" exact component={User} />
        <Route path="/" render={ ()=> { return main()} } />
      </Switch>
    );

  return (
    <div className={classes.App}>
      <Header 
        logout={handleLogout} 
        login={handleLogin}
        handleUserAccountInfo={handleUserAccountInfo}
      />
      <div>
        {content}
      </div>
      <Footer/>
    </div>
  )
  
}

// need use withLD here so the LD component can receive the user property !
export default withRouter(withLD(Layout));