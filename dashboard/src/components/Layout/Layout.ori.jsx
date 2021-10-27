
import React, { Component } from "react";

import { Grid } from "@material-ui/core";
import CircularProgress from '@material-ui/core/CircularProgress';
import { Route, withRouter, Switch } from "react-router-dom";

import {connect} from 'react-redux';



// import ContentPanel from "../Contents/ContentPanel";
import Box from "@material-ui/core/Box";
import Search from "../Contents/Search/Search";
// import SystemInfo from "../Contents/SystemInfo/SystemInfoDemo";
import SystemInfo from "../Contents/SystemInfo/SystemInfo";
import JobInfo from "../Contents/JobInfo/JobInfo";
import Workflow from "../Contents/Workflow/Workflow";
import User from "../Contents/User/User";
import Advanced from "../Contents/PipelineBuilder/PipelineBuilder";
import MyJobs from "../Contents/JobInfo/MyJobs";

import BlastDB from "../Contents/BlastDB/BlastDB";

import Header from "../Header/Header";
import Footer from "../Footer/Footer";
import * as actions from  "../../store/actions";

import Controls from "../Menu/Controls";
import withLD from "../hoc/withLD";

import {cookieGet, cookieSet, setSSOReturnURL, getSSOHash, CK_AVATAR_ID} from "../../utils/cookie";
import { ssoURL } from "../../utils/sso";

import classes from "./Layout.module.css";

class Layout extends Component {
  pageMap = [
    {label: "My Jobs", path:"/myjobs", component: MyJobs},
    {label: "Site Search", path: "/", component: Search}, 
    {label: "System info", path: "/sysinfo", component: SystemInfo}, 
    {label: "Job info", path: "/jobinfo", component: JobInfo}, 
    {label: "Workflow Catalog", path: "/workflow", component: Workflow},
    {label: "BlastDB Catalog", path: "/blastdb", component: BlastDB},
    {label: "Pipeline Builder", path:"/builder", component:Advanced}
  ]

  componentDidMount() {
    this.updateAvatarId();

    const sso_hash = getSSOHash();
    // console.error(sso_session, "[Layout]");
    
    if (sso_hash){
      this.props.onAuth(sso_hash);
    } 
    // this.props.history.push(this.props.location.pathname);  // keep the current URL
  }

  //sync avarta id
  updateAvatarId(){
    let aid = cookieGet(CK_AVATAR_ID);
    if (!aid) {
      cookieSet(CK_AVATAR_ID, this.props.avatarId);
    } else {
      aid = parseInt(aid);
      this.props.updateAvatarId(aid);
    }
  }

  handleLogout = () => {
    setSSOReturnURL();
    window.location.href = `${ssoURL}/destroy`;
    this.props.onLogout();
  }

  main = () => {
    const pmap = this.pageMap.map( p => ({...p}) ); //make a copy ...

    // console.log(pmap, "[routes] map")
    // need to pass flags in props down to page component ...
    const routes = pmap.map( p => <Route key={p.label} path={p.path} exact 
      render={ () => { return <p.component {...this.props}/>} } />)

  
    // console.log(routes, "[routes] comp")
    return (
      <>
        <Grid container justify="flex-start">
          <Grid item xl={2}>
            <Controls controls={pmap.map(p=>[p.label, p.path])} {...this.props}/>
          </Grid>  
          <Grid item xs={10} sm={10} md={9} lg={9} xl={9}>
            <div>
              <Box p={2}>
                <Switch>
                  { routes }
                </Switch>
              </Box>
            </div>
          </Grid>
        </Grid>
      </>
    )
  }


  render() { 
 
    const content = this.props.loading ? <CircularProgress /> :
      (
        <Switch>
          <Route path="/user" exact component={User} />
          <Route path="/" render={ ()=> { return this.main()} } />
        </Switch>
      );

    return (
      <div className={classes.App}>
        <Header logout={this.handleLogout} />
        <div>
          {content}
        </div>
        <Footer/>
      </div>
    )
  }
}

const mapStateToProps = state => {
  return {
    loading: state.auth.loading,
    avatarId : state.main.avatarId,
    user: state.auth.user
  }
}

const mapDispatchToProps = dispatch => {
  return {
    updateAvatarId: (aid) => dispatch(actions.updateUserAvatarId(aid)),
    onAuth: (sso_token) => dispatch(actions.auth(sso_token)),
    onLogout: () => dispatch(actions.logout())
  }
}

// need use withLD here so the LD component can receive the user property !
export default connect(mapStateToProps, mapDispatchToProps)(withRouter(withLD(Layout)));