import React from "react";

import AppBar from "@material-ui/core/AppBar";
import Toolbar from "@material-ui/core/Toolbar";
import Typography from "@material-ui/core/Typography";
import Button from "@material-ui/core/Button";
import IconButton from "@material-ui/core/IconButton";

import AccountCircleIcon from '@material-ui/icons/AccountCircle';

import { Alert, Skeleton } from "@material-ui/lab";

import { Link } from "react-router-dom";

import { connect } from "react-redux";

import classes from "./Header.module.css";
import logo from "../../assets/images/jaws.png";
// import { avatars } from "../../utils/avatars";
import UserAvatar from "./UserAvatar";

// import { handleLogin } from "../../utils/sso";

const Header = ({
  error, 
  loading, 
  user, 
  login, 
  logout, 
  avatarId,
  handleUserAccountInfo,
}) => {

    let userui = ( <Button onClick={login}>
                      <AccountCircleIcon className={classes.defaultUserIcon}/>
                    </Button> )
    if (error){
        userui = ( <Alert
                    severity="error"
                    variant="filled"
                    size="10px"
                    className={classes.Alert}
                    onClose={()=>{ logout() } }
                    >User info error</Alert>
                )
    } else if (loading){
      userui = <Skeleton animation="wave" style={{ width: "200px", height: "24px" }} />
    } else if (user) {
      const uobj = user;
      userui =  ( <ul className={classes.Ul}>
                    <li className={classes.Li} title={`${uobj.first_name} ${uobj.last_name}`}>
                      <UserAvatar
                        user={uobj}
                        handleLogout={logout}
                        avatarId={avatarId}
                        handleUserAccountInfo={handleUserAccountInfo}
                      />
                    </li>
                  </ul>
                )
    }

    return (     
      <div className={classes.Root}>
        <AppBar position="static" color="transparent">
          <Toolbar>
            <IconButton edge="start" className={classes.MenuButton} color="inherit" aria-label="menu">
              <Link to="/">
                <img src={logo} className={classes.Logo} alt="JAWS" />
              </Link>
            </IconButton>
            <Typography variant="h6" className={classes.Title} id="myDIV">
              <Link to="/"
                  className={classes.NoDeco}>
                JGI JAWs
              </Link>
              <Link to="//jaws-docs.readthedocs.io" 
                  className={classes.headerLink} 
                  target="_blank" 
                  rel="external noopener noreferrer">
                Docs
              </Link>
              <Link to="//grafana.jgi.lbl.gov/d/S7d2rMPMk/jaws-status?orgId=1&refresh=1m" 
                  className={classes.headerLink} 
                  target="_blank" 
                  rel="external noopener noreferrer">
                JAWS Services
              </Link>
            </Typography>
            {userui}
          </Toolbar>
        </AppBar>
      </div> 
    );
}
 
const mapStateToProps = state => {
  return {
    avatarId: state.main.avatarId,
    user: state.auth.user,
    loading: state.auth.loading,
    error: state.auth.error
  }
}
export default connect(mapStateToProps)(Header);

