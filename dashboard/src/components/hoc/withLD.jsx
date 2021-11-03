import React, { Component } from "react";
import { connect } from "react-redux";

// Ref : https://docs.launchdarkly.com/sdk/client-side/react

import { withLDProvider, withLDConsumer } from "launchdarkly-react-client-sdk";
// import { lightBlue } from "@material-ui/core/colors";

const withLaunchdarkly = WrappedComponent => {
  class LaunchDarklyHOC extends Component {
    componentDidUpdate = prevProps => {
      let pemail = null; // user in prevProps
      let cemail = null; // user in current props

      // console.log(this.props, "[withLD] didUpdate");

      if (prevProps.user?.email_address){
        pemail = prevProps.user.email_address;
      }

      if(this.props.user?.email_address){
        cemail = this.props.user.email_address;
      }
      // console.log(pemail, cemail, this.props, "[withLD] hoc users");
      const { ldClient } = this.props;
      if (pemail !== cemail && ldClient && this.props.user) {
        // console.log(this.props, "[withLD] this.prop")
        // User changed, update current LD user
        const user = {
          key: cemail,
          name: `${this.props.user.first_name} ${this.props.user.last_name}`
        }
        // console.log(user, "[withLD] user")
        ldClient.identify(user);
      }
    };

    render() {
      // console.log(this.props, "[withLD] render");
      return <WrappedComponent {...this.props} />;
    }
  }

  const mapStateToProps = state => {
    return {
      user: state.auth,
    }
  }

  // JDP LD account, jaws-dashboard project :
  // React SDK config 
  const ldClientKey = process.env.NODE_ENV === "production"
          ? process.env.REACT_APP_LD_CLIENT_KEY_PRODUCTION
          : process.env.REACT_APP_LD_CLIENT_KEY_TEST;
  const lDConfig = { clientSideID: ldClientKey };

  return withLDProvider(lDConfig)(withLDConsumer()(connect(mapStateToProps)(LaunchDarklyHOC)));
};



export default withLaunchdarkly;
